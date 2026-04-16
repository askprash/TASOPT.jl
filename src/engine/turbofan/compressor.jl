"""
Compressor component type for a TASOPT turbofan.

Models an axial compressor stage (fan, LPC, or HPC) from its inlet station to
its outlet station.  The composed map subtype is the gridded-interpolation
`CompressorMap` from `map_functions.jl`, which stores 2D tables for corrected
mass flow, pressure ratio, and polytropic efficiency over a (speed, R-line)
grid.
"""

# ---------------------------------------------------------------------------
# Compressor struct  —  typed component with composed CompressorMap subtype
# ---------------------------------------------------------------------------

"""
    Compressor{T<:AbstractFloat}

Parametric compressor component holding design-point anchors, an efficiency
floor, and a reference to a composed `CompressorMap` subtype.

Instantiate one `Compressor` per spool stage:

```julia
fan = Compressor(pifD,  mbfD,  NbfD,  epf0,  T(0.60), FanMap)
lpc = Compressor(pilcD, mblcD, NblcD, eplc0, T(0.70), LPCMap)
hpc = Compressor(pihcD, mbhcD, NbhcD, ephc0, T(0.70), HPCMap)
```

## Fields

| Field       | Unit  | Description                                                    |
|:-----------|:------|:---------------------------------------------------------------|
| `piD`       | —     | Design compression pressure ratio  pt_out / pt_in  (> 1)      |
| `mbD`       | kg/s  | Design corrected mass flow  (= ṁ√(T/Tref) / (p/pref))         |
| `NbD`       | —     | Design corrected spool speed  (= N / √(T/Tref))               |
| `epol0`     | —     | Maximum polytropic efficiency at the design point  (< 1)      |
| `epol_min`  | —     | Efficiency floor; output is clamped to this minimum           |
| `map`       | —     | Composed `CompressorMap` holding 2D interpolation tables       |
| `Ng`        | —     | Warm-start speed guess for the map-inversion solver            |
| `Rg`        | —     | Warm-start R-line guess for the map-inversion solver           |

## Station convention

This type is agnostic of station numbers; the caller (Newton driver or
`tfoper!`) is responsible for supplying the correct `FlowStation` arguments.

| Stage | Inlet station | Outlet station | Label                     |
|:------|:-------------|:--------------|:--------------------------|
| Fan   | 2            | 21             | fan compressor            |
| LPC   | 19c          | 25             | low-pressure compressor   |
| HPC   | 25c          | 3              | high-pressure compressor  |
"""
mutable struct Compressor{T<:AbstractFloat}
    piD      ::T             # design compression pressure ratio  pt_out/pt_in   [—]
    mbD      ::T             # design corrected mass flow                         [kg/s]
    NbD      ::T             # design corrected spool speed                       [—]
    epol0    ::T             # max polytropic efficiency                          [—]
    epol_min ::T             # efficiency floor  (≥ 0, < 1)                      [—]
    map      ::CompressorMap # composed map subtype  (FanMap, LPCMap, or HPCMap)
    Ng       ::Float64       # map-solver speed guess    (mutable solver state)  [—]
    Rg       ::Float64       # map-solver R-line guess   (mutable solver state)  [—]
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    Compressor(piD, mbD, NbD, epol0, epol_min, map; Ng=0.5, Rg=2.0)

Convenience constructor for a `Compressor` component.

`piD`, `mbD`, `NbD`, `epol0`, and `epol_min` set the design anchors and the
efficiency floor; `map` is the `CompressorMap` instance (one of `FanMap`,
`LPCMap`, `HPCMap` from `maps.jl`).  `Ng` and `Rg` are optional warm-start
guesses for the internal NLsolve map-inversion call (default: 0.5 and 2.0,
near the map centre).
"""
function Compressor(
    piD      ::T,
    mbD      ::T,
    NbD      ::T,
    epol0    ::T,
    epol_min ::T,
    map      ::CompressorMap;
    Ng       ::Float64 = 0.5,
    Rg       ::Float64 = 2.0,
) where {T<:AbstractFloat}
    Compressor{T}(piD, mbD, NbD, epol0, epol_min, map, Ng, Rg)
end

# ---------------------------------------------------------------------------
# compressor_efficiency  —  corrected speed and polytropic efficiency from map
# ---------------------------------------------------------------------------

"""
    compressor_efficiency(comp, pi, mb)
      -> (Nb, epol, Nb_pi, Nb_mb, epol_pi, epol_mb)

Return the corrected spool speed, polytropic efficiency, and their partial
derivatives for the compressor stage `comp` at the off-design operating point
`(pi, mb)`.

This is a thin wrapper around `calculate_compressor_speed_and_efficiency`,
binding the design-point anchors from `comp` and applying two post-processing
corrections:

1. **Efficiency floor**: if the map returns an efficiency below `comp.epol_min`,
   the efficiency is clamped to the floor and both derivatives are zeroed.
2. **Windmilling inversion**: if `pi < 1` (pressure ratio below unity), the
   compressor is operating in reverse (turbine-like windmilling mode).
   Following the TASOPT convention, the efficiency is replaced by its
   reciprocal and the derivatives are adjusted via the chain rule:
   `d(1/epol)/dpi = −(1/epol²) · depol/dpi`.

The map-solver warm-start guesses `comp.Ng` and `comp.Rg` are updated in-place
after each call so that subsequent evaluations at nearby operating points
converge faster.

## Arguments

| Argument | Unit  | Description                                           |
|:---------|:------|:------------------------------------------------------|
| `comp`   | —     | `Compressor` component                                |
| `pi`     | —     | Off-design compression pressure ratio  pt_out / pt_in |
| `mb`     | kg/s  | Off-design corrected mass flow                        |

## Returns

A 6-tuple `(Nb, epol, Nb_pi, Nb_mb, epol_pi, epol_mb)`:

| Value      | Description                                       |
|:-----------|:--------------------------------------------------|
| `Nb`       | Corrected spool speed                             |
| `epol`     | Polytropic efficiency (after floor + windmilling) |
| `Nb_pi`    | ∂Nb / ∂pi                                         |
| `Nb_mb`    | ∂Nb / ∂mb                                         |
| `epol_pi`  | ∂epol / ∂pi  (after floor + windmilling)          |
| `epol_mb`  | ∂epol / ∂mb  (after floor + windmilling)          |

## Design-point identity

At the design operating point (`pi = piD`, `mb = mbD`), the map interpolates
at its reference conditions and the scaled efficiency equals `epol0`.
"""
function compressor_efficiency(
    comp ::Compressor{T},
    pi   ::T,
    mb   ::T,
) where {T<:AbstractFloat}
    # Map-inversion solver requires Float64; convert T if needed
    Nb, epol, Nb_pi, Nb_mb, epol_pi, epol_mb, N, R =
        calculate_compressor_speed_and_efficiency(
            comp.map,
            Float64(pi), Float64(mb),
            Float64(comp.piD), Float64(comp.mbD), Float64(comp.NbD), Float64(comp.epol0),
            Ng = comp.Ng, Rg = comp.Rg,
        )

    # Update warm-start hints for the next Newton iteration
    comp.Ng = N
    comp.Rg = R

    epol_floor = Float64(comp.epol_min)

    # 1. Efficiency floor: clamp if below minimum
    if epol < epol_floor
        epol    = epol_floor
        epol_pi = 0.0
        epol_mb = 0.0
    end

    # 2. Windmilling inversion (pi < 1 → compressor acting as turbine)
    #    TASOPT convention: use reciprocal efficiency in gas_pratd
    if Float64(pi) < 1.0
        epol_pi = (-1.0 / epol^2) * epol_pi
        epol_mb = (-1.0 / epol^2) * epol_mb
        epol    = 1.0 / epol
    end

    return T(Nb), T(epol), T(Nb_pi), T(Nb_mb), T(epol_pi), T(epol_mb)
end

# ---------------------------------------------------------------------------
# compressor_exit!  —  compute outlet FlowStation from inlet + compression
# ---------------------------------------------------------------------------

"""
    compressor_exit!(st_out, st_in, pi, epol) -> st_out

Compute the outlet `FlowStation` of a compressor stage given the inlet state,
the compression pressure ratio, and the stage polytropic efficiency.

The outlet total state is obtained by calling `gas_pratd` with the pressure
ratio `pi` and efficiency `epol`.  Only the six primary thermodynamic scalars
(`pt, Tt, ht, st, cpt, Rt`) and the gas composition vector `alpha` are written
to `st_out`; static-state fields (`Ts, ps, u, …`) are **not** set here because
they depend on the local Mach number, which is determined by downstream flow
matching.

## Arguments

| Argument | Unit  | Description                                        |
|:---------|:------|:---------------------------------------------------|
| `st_out` | —     | Outlet `FlowStation` (mutated in place)            |
| `st_in`  | —     | Inlet `FlowStation` (read-only)                    |
| `pi`     | —     | Compression pressure ratio  pt_out / pt_in  (> 1)  |
| `epol`   | —     | Stage polytropic efficiency  (0 < epol < 1)        |

## Physics

`gas_pratd` solves iteratively for the outlet temperature `Tt_out` consistent
with the polytropic process:

    (s_out − s_in) / R = log(pi) / epol

and computes the outlet total pressure as `pt_out = pt_in · pi`.  For an ideal
compressor `epol → 1`, the process approaches isentropic; for `epol < 1`,
additional entropy is generated and the total-temperature rise is larger than
the isentropic value for the same pressure ratio (less efficient compression).

## Returns

`st_out` (mutated in place).
"""
function compressor_exit!(
    st_out ::FlowStation{T},
    st_in  ::FlowStation{T},
    pi     ::T,
    epol   ::T,
) where {T<:AbstractFloat}
    # gas_pratd returns (pt, Tt, ht, st, cpt, Rt, pt_pt_in, ...derivatives...)
    res = gas_pratd(
        st_in.alpha, 5,
        st_in.pt, st_in.Tt, st_in.ht, st_in.st, st_in.cpt, st_in.Rt,
        pi, epol,
    )

    st_out.pt    = res[1]
    st_out.Tt    = res[2]
    st_out.ht    = res[3]
    st_out.st    = res[4]
    st_out.cpt   = res[5]
    st_out.Rt    = res[6]
    st_out.alpha = st_in.alpha   # composition unchanged through compressor

    return st_out
end
