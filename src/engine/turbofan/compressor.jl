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
| `windmilling` | —   | Enable windmilling inversion when `pi < 1` (fan only)         |
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
    piD        ::T             # design compression pressure ratio  pt_out/pt_in   [—]
    mbD        ::T             # design corrected mass flow                         [kg/s]
    NbD        ::T             # design corrected spool speed                       [—]
    epol0      ::T             # max polytropic efficiency                          [—]
    epol_min   ::T             # efficiency floor  (≥ 0, < 1)                      [—]
    map        ::CompressorMap # composed map subtype  (FanMap, LPCMap, or HPCMap)
    windmilling::Bool          # enable pi<1 efficiency inversion (fan only)       [—]
    Ng         ::Float64       # map-solver speed guess    (mutable solver state)  [—]
    Rg         ::Float64       # map-solver R-line guess   (mutable solver state)  [—]
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    Compressor(piD, mbD, NbD, epol0, epol_min, map; windmilling=false, Ng=0.5, Rg=2.0)

Convenience constructor for a `Compressor` component.

`piD`, `mbD`, `NbD`, `epol0`, and `epol_min` set the design anchors and the
efficiency floor; `map` is the `CompressorMap` instance (one of `FanMap`,
`LPCMap`, `HPCMap` from `maps.jl`).  Set `windmilling=true` for the fan stage
to enable the `pi < 1` efficiency inversion (upstream TASOPT applies this only
to the fan).  `Ng` and `Rg` are optional warm-start guesses for the internal
NLsolve map-inversion call (default: 0.5 and 2.0, near the map centre).
"""
function Compressor(
    piD      ::T,
    mbD      ::T,
    NbD      ::T,
    epol0    ::T,
    epol_min ::T,
    map      ::CompressorMap;
    windmilling ::Bool   = false,
    Ng       ::Float64 = 0.5,
    Rg       ::Float64 = 2.0,
) where {T<:AbstractFloat}
    Compressor{T}(piD, mbD, NbD, epol0, epol_min, map, windmilling, Ng, Rg)
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
    #    TASOPT convention: use reciprocal efficiency in gas_pratd.
    #    Upstream applies this only to the fan; LPC/HPC never windmill.
    if comp.windmilling && Float64(pi) < 1.0
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
# ---------------------------------------------------------------------------
# compressor_Nb_residual  —  corrected-speed map-match residual
# ---------------------------------------------------------------------------

"""
    compressor_Nb_residual(comp, pi, mb, Nb_target) -> (r, Nb_pi, Nb_mb)

Map-match residual for the spool-speed coupling constraint.

In a TASOPT turbofan the fan and LPC share a shaft (via gear ratio `Gearf`).
At a converged off-design operating point the corrected spool speeds must
satisfy:

    Gearf × √(Tt2/Tref) × Nb_fan = √(Tt19c/Tref) × Nb_lpc

The per-compressor primitive is:

    r = Nb_map(pi, mb) − Nb_target

where `Nb_target` is the corrected speed imposed by the complementary spool
member (e.g. `Nl × trl / (Gearf × trf)` for the fan side).  The Newton
driver assembles the full row-1 residual from two calls:

    res[1] = Gearf × trf × r_fan                # r_fan with Nb_target = Nl×trl/(Gearf×trf)
           = Gearf × trf × Nb_fan − trl × Nb_lpc

At any converged Newton point the shaft coupling is satisfied and `r = 0`.

The partial derivatives of r with respect to `pi` and `mb` are those of
`Nb_map` (from `compressor_efficiency`); the derivative with respect to
`Nb_target` is always −1.

## Arguments

| Argument    | Unit  | Description                                            |
|:-----------|:------|:-------------------------------------------------------|
| `comp`      | —     | `Compressor` component                                 |
| `pi`        | —     | Off-design compression pressure ratio  pt_out / pt_in  |
| `mb`        | kg/s  | Off-design corrected mass flow                         |
| `Nb_target` | —     | Target corrected spool speed from the shaft coupling   |

## Returns

A 3-tuple `(r, Nb_pi, Nb_mb)`:

| Value    | Description                              |
|:---------|:-----------------------------------------|
| `r`      | `Nb_map(pi, mb) − Nb_target`             |
| `Nb_pi`  | ∂Nb_map / ∂pi                            |
| `Nb_mb`  | ∂Nb_map / ∂mb                            |

## Design-point identity

At the design operating point (`pi = piD`, `mb = mbD`) the map scaling
returns `Nb_map = NbD` exactly (see `calculate_compressor_speed_and_efficiency`).
Therefore, calling with `Nb_target = comp.NbD` gives `r = 0` exactly.
"""
function compressor_Nb_residual(
    comp      ::Compressor{T},
    pi        ::T,
    mb        ::T,
    Nb_target ::T,
) where {T<:AbstractFloat}
    Nb, _, Nb_pi, Nb_mb, _, _ = compressor_efficiency(comp, pi, mb)
    r = Nb - Nb_target
    return r, Nb_pi, Nb_mb
end

# ---------------------------------------------------------------------------
# compressor_exit!  —  compute outlet FlowStation from inlet + compression
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# compressor_pratd  —  outlet state + analytic Jacobian block
# ---------------------------------------------------------------------------

"""
    compressor_pratd(comp, alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, pi, mb)
      -> (pt_out, Tt_out, ht_out, st_out, cpt_out, Rt_out,
          pt_out_pt_in,
          pt_out_st_in, Tt_out_st_in, ht_out_st_in, st_out_st_in,
          pt_out_pi, Tt_out_pi, ht_out_pi, st_out_pi,
          pt_out_mb, Tt_out_mb, ht_out_mb, st_out_mb,
          Nb, Nb_pi, Nb_mb, epol, epol_pi, epol_mb)

Compute the outlet total state of a compressor stage and its analytic partial
derivatives for use in the Newton Jacobian.

This function combines `compressor_efficiency` (map inversion) with `gas_pratd`
(isentropic-in-efficiency compression) and applies the efficiency chain rule to
produce the full derivative set needed by the Newton driver.

## Arguments

| Argument  | Unit   | Description                                         |
|:----------|:-------|:----------------------------------------------------|
| `comp`    | —      | `Compressor` component (design anchors + map)       |
| `alpha`   | —      | Gas composition mass fractions                      |
| `nair`    | —      | Number of gas species                               |
| `pt_in`   | Pa     | Inlet total pressure                                |
| `Tt_in`   | K      | Inlet total temperature                             |
| `ht_in`   | J/kg   | Inlet total enthalpy                                |
| `st_in`   | J/kg·K | Inlet entropy complement                            |
| `cpt_in`  | J/kg·K | Inlet specific heat                                 |
| `Rt_in`   | J/kg·K | Inlet gas constant                                  |
| `pi`      | —      | Compression ratio  pt_out / pt_in  (Newton variable)|
| `mb`      | kg/s   | Corrected mass flow  (Newton variable)              |

## Returns

A 26-tuple. Derivative naming convention: `X_Y` means ∂X/∂Y.

| Value               | Description                                        |
|:--------------------|:---------------------------------------------------|
| `pt_out … Rt_out`   | Outlet total state (6 scalars)                     |
| `pt_out_pt_in`      | ∂pt_out/∂pt_in  (= `pi`)                           |
| `pt_out_st_in`      | ∂pt_out/∂st_in  (= 0; pressure is pt_in·pi)        |
| `Tt_out_st_in …`    | ∂Tt_out/∂st_in, ∂ht_out/∂st_in, ∂st_out/∂st_in   |
| `pt_out_pi …`       | ∂X/∂pi with efficiency chain rule applied (4)      |
| `pt_out_mb …`       | ∂X/∂mb, efficiency-chain only; add `pt_out_pt_in · ∂pt_in/∂mb` in caller (4) |
| `Nb, Nb_pi, Nb_mb`  | Corrected spool speed and its derivatives (shaft Jacobian) |
| `epol, epol_pi, epol_mb` | Polytropic efficiency and its derivatives     |

## Notes

The `pt_out_mb` block captures only the efficiency-chain contribution
`∂pt_out/∂epol · ∂epol/∂mb`.  The caller must add the upstream total-pressure
chain `pt_out_pt_in · ∂pt_in/∂mb` to form the complete `∂pt_out/∂mb`.
This separation lets the function remain agnostic of which Newton variable
drives `pt_in` (e.g. inlet Mach number or corrected mass flow of an upstream
compressor).
"""
function compressor_pratd(
    comp   ::Compressor{T},
    alpha, nair,
    pt_in  ::T, Tt_in ::T, ht_in ::T, st_in ::T, cpt_in ::T, Rt_in ::T,
    pi     ::T,
    mb     ::T,
) where {T<:AbstractFloat}

    # 1. Map inversion: corrected speed + polytropic efficiency and derivatives
    Nb, epol, Nb_pi, Nb_mb, epol_pi, epol_mb = compressor_efficiency(comp, pi, mb)

    # 2. Isentropic-in-efficiency compression (gas_pratd requires Float64)
    res = gas_pratd(alpha, nair,
        Float64(pt_in), Float64(Tt_in), Float64(ht_in), Float64(st_in),
        Float64(cpt_in), Float64(Rt_in), Float64(pi), Float64(epol))

    pt_out  = T(res[1]);  Tt_out  = T(res[2])
    ht_out  = T(res[3]);  st_out  = T(res[4])
    cpt_out = T(res[5]);  Rt_out  = T(res[6])

    pt_out_pt_in = T(res[7])   # = pi  (pt_out = pt_in · pi)
    # res[8] = p_so = 0 always (pressure doesn't depend on inlet entropy)
    pt_out_st_in  = T(res[8])  # 0
    Tt_out_st_in  = T(res[9])
    ht_out_st_in  = T(res[10])
    st_out_st_in  = T(res[11])

    pt_out_pi_raw = T(res[12]);  Tt_out_pi_raw = T(res[13])
    ht_out_pi_raw = T(res[14]);  st_out_pi_raw = T(res[15])

    pt_out_epol = T(res[16]);  Tt_out_epol = T(res[17])
    ht_out_epol = T(res[18]);  st_out_epol = T(res[19])

    # 3. Chain rule: apply efficiency dependence on pi
    pt_out_pi = pt_out_epol * epol_pi + pt_out_pi_raw
    Tt_out_pi = Tt_out_epol * epol_pi + Tt_out_pi_raw
    ht_out_pi = ht_out_epol * epol_pi + ht_out_pi_raw
    st_out_pi = st_out_epol * epol_pi + st_out_pi_raw

    # 4. mb direction: efficiency chain only
    #    Caller adds pt_out_pt_in * (∂pt_in/∂mb) for the upstream pressure contribution
    pt_out_mb = pt_out_epol * epol_mb
    Tt_out_mb = Tt_out_epol * epol_mb
    ht_out_mb = ht_out_epol * epol_mb
    st_out_mb = st_out_epol * epol_mb

    return pt_out, Tt_out, ht_out, st_out, cpt_out, Rt_out,
           pt_out_pt_in,
           pt_out_st_in, Tt_out_st_in, ht_out_st_in, st_out_st_in,
           pt_out_pi, Tt_out_pi, ht_out_pi, st_out_pi,
           pt_out_mb, Tt_out_mb, ht_out_mb, st_out_mb,
           Nb, Nb_pi, Nb_mb, epol, epol_pi, epol_mb
end
