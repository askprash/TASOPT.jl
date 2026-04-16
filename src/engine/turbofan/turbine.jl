"""
Turbine component type for a TASOPT turbofan.

Models an axial turbine stage (HPT or LPT) from its inlet station to its
outlet station.  The map subtype is the standard quadratic polynomial fit
from `etmap` (tfmap.jl), parameterised by two sensitivity constants `pcon`
and `Ncon` stored in the `Tmap` field.
"""

# ---------------------------------------------------------------------------
# TurbineMap struct  —  composed map subtype
# ---------------------------------------------------------------------------

"""
    TurbineMap{T<:AbstractFloat}

Quadratic penalty map for a turbine stage, encoding how polytropic efficiency
degrades away from the design point.

The efficiency model (Drela TASOPT tfmap.jl, function `etmap`) is:

    prat  = (T_in / T_out)^(cpt / (Rt · ep0))   # isentropic-equivalent expansion ratio
    Nmb   = Nb · mb                              # speed–flow product
    NmbD  = NbD · mbD                            # design-point speed–flow product

    ept = ep0 · (1 − pcon·(1 − prat/piD)² − Ncon·(1 − Nmb/NmbD)²)

## Fields

| Field  | Unit | Description                                              |
|:-------|:-----|:---------------------------------------------------------|
| `pcon` | —    | Pressure-sensitivity constant; controls efficiency droop |
|        |      | as the expansion ratio departs from its design value     |
| `Ncon` | —    | Speed-sensitivity constant; controls efficiency droop    |
|        |      | as the corrected speed–flow product departs from design  |
"""
struct TurbineMap{T<:AbstractFloat}
    pcon::T   # pressure-sensitivity constant  [—]
    Ncon::T   # speed-sensitivity constant     [—]
end

"""
    TurbineMap(v::AbstractVector)

Construct a `TurbineMap` from a 2-element vector `[pcon, Ncon]`, matching the
layout of the legacy `Tmaph`/`Tmapl` constants in `maps.jl`.
"""
function TurbineMap(v::AbstractVector{T}) where {T<:AbstractFloat}
    @assert length(v) == 2 "TurbineMap vector must have exactly 2 elements [pcon, Ncon]"
    TurbineMap{T}(v[1], v[2])
end

# ---------------------------------------------------------------------------
# Turbine struct
# ---------------------------------------------------------------------------

"""
    Turbine{T<:AbstractFloat}

Parametric turbine component holding design-point anchors and a composed
`TurbineMap` subtype.

Instantiate one `Turbine` per spool stage:

```julia
hpt = Turbine(pihtD, mbhtD, NbhtD, epht0; map=TurbineMap(Tmaph))
lpt = Turbine(piltD, mbltD, NbltD, eplt0; map=TurbineMap(Tmapl))
```

## Fields

| Field  | Unit          | Description                                               |
|:-------|:--------------|:----------------------------------------------------------|
| `piD`  | —             | Design expansion pressure ratio  pt_in / pt_out (> 1)    |
| `mbD`  | kg/s          | Design corrected mass flow  (= ṁ√(T/Tref) / (p/pref) )   |
| `NbD`  | —             | Design corrected spool speed  (= N / √(T/Tref) )         |
| `ep0`  | —             | Maximum polytropic efficiency at the design point (< 1)  |
| `map`  | —             | Composed `TurbineMap{T}` holding the penalty constants    |

## Station convention

This type is agnostic of station numbers; the caller (Newton driver or
`tfoper!`) is responsible for supplying the correct `FlowStation` arguments.

| Stage | Inlet station | Outlet station | Typical label |
|:------|:-------------|:--------------|:--------------|
| HPT   | 41 (post-combustor + cooling mixing) | 45 | high-pressure turbine |
| LPT   | 45           | 49             | low-pressure turbine  |
"""
mutable struct Turbine{T<:AbstractFloat}
    piD ::T              # design expansion pressure ratio  pt_in/pt_out  [—]
    mbD ::T              # design corrected mass flow                      [kg/s]
    NbD ::T              # design corrected spool speed                    [—]
    ep0 ::T              # max polytropic efficiency                       [—]
    map ::TurbineMap{T}  # composed map subtype
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    Turbine(piD, mbD, NbD, ep0; map=TurbineMap(0.15, 0.15))

Convenience constructor for a `Turbine` component.  `piD`, `mbD`, `NbD`, and
`ep0` are the design-point anchors; `map` defaults to the standard TASOPT
turbine map constants (pcon = Ncon = 0.15).
"""
function Turbine(
    piD ::T,
    mbD ::T,
    NbD ::T,
    ep0 ::T;
    map ::TurbineMap{T} = TurbineMap{T}(T(0.15), T(0.15)),
) where {T<:AbstractFloat}
    Turbine{T}(piD, mbD, NbD, ep0, map)
end

# ---------------------------------------------------------------------------
# turbine_efficiency  —  polytropic efficiency from the composed map
# ---------------------------------------------------------------------------

"""
    turbine_efficiency(turb, dh, mb, Nb, Tt, cpt, Rt)
      -> (ept, ept_dh, ept_mb, ept_Nb, ept_Tt, ept_cpt, ept_Rt)

Return the polytropic efficiency and its analytic partial derivatives for the
turbine stage `turb`, given the off-design operating point.

This is a thin wrapper around the low-level `etmap` function from `tfmap.jl`,
binding the design-point anchors and map constants from `turb`.

## Arguments

| Argument | Unit  | Description                                            |
|:---------|:------|:-------------------------------------------------------|
| `turb`   | —     | `Turbine` component                                    |
| `dh`     | J/kg  | Specific enthalpy drop across the stage (< 0)         |
| `mb`     | kg/s  | Off-design corrected mass flow                         |
| `Nb`     | —     | Off-design corrected spool speed                       |
| `Tt`     | K     | Inlet total temperature                                |
| `cpt`    | J/kg·K| Inlet total specific heat                              |
| `Rt`     | J/kg·K| Inlet gas constant                                     |

## Returns

A 7-tuple `(ept, ept_dh, ept_mb, ept_Nb, ept_Tt, ept_cpt, ept_Rt)` matching
the return convention of `etmap`.  The analytic derivatives are needed by the
Newton driver (tasopt-j9l.31) for Jacobian assembly.

## Design-point identity

At the design operating point (`mb = mbD`, `Nb = NbD`, `dh` matches `piD`),
the efficiency equals `ep0` because both penalty terms vanish:

    prat/piD = 1   →   pcon-term = 0
    Nmb/NmbD = 1   →   Ncon-term = 0
    ept = ep0 · (1 − 0 − 0) = ep0
"""
function turbine_efficiency(
    turb ::Turbine{T},
    dh   ::T,
    mb   ::T,
    Nb   ::T,
    Tt   ::T,
    cpt  ::T,
    Rt   ::T,
) where {T<:AbstractFloat}
    etmap(dh, mb, Nb, turb.piD, turb.mbD, turb.NbD, turb.ep0,
          (turb.map.pcon, turb.map.Ncon), Tt, cpt, Rt)
end

# ---------------------------------------------------------------------------
# turbine_exit!  —  compute outlet FlowStation from inlet + work extraction
# ---------------------------------------------------------------------------

"""
    turbine_exit!(st_out, st_in, dh, ept) -> st_out

Compute the outlet `FlowStation` of a turbine stage given the inlet state, a
specific enthalpy drop, and the stage polytropic efficiency.

The outlet total state is obtained by calling `gas_delhd` with the inverse
efficiency `epi = 1/ept` as the polytropic exponent.  Only the six primary
thermodynamic scalars (pt, Tt, ht, st, cpt, Rt) and the gas composition
vector `alpha` are written to `st_out`; static-state fields (Ts, ps, u, …)
are **not** set here because they depend on the local Mach number, which is
determined by downstream flow matching.

## Arguments

| Argument | Unit  | Description                                       |
|:---------|:------|:--------------------------------------------------|
| `st_out` | —     | Outlet `FlowStation` (mutated in place)           |
| `st_in`  | —     | Inlet `FlowStation` (read-only)                   |
| `dh`     | J/kg  | Specific enthalpy change (< 0 for work extraction)|
| `ept`    | —     | Stage polytropic efficiency (0 < ept < 1)         |

## Physics

`gas_delhd` solves iteratively for the outlet temperature `Tt_out` consistent
with `h(alpha, Tt_out) = h_in + dh`, then computes the outlet total pressure
using the polytropic relation:

    pt_out = pt_in · exp(epi · (s_out − s_in) / Rt_in)

where `epi = 1/ept`.  For an ideal turbine `ept → 1`, the process
approaches isentropic; for `ept < 1`, additional entropy is generated and
the total-pressure drop is smaller than the isentropic value for the same
work extraction (i.e., less efficient expansion).

## Returns

`st_out` (mutated in place).
"""
# ---------------------------------------------------------------------------
# turbine_mb_residual  —  map-match residual for vertical-line turbine matching
# ---------------------------------------------------------------------------

"""
    turbine_mb_residual(turb, mb) -> (r, dr_dmb)

Map-match residual for the vertical-line turbine constraint.

Turbines in TASOPT are matched on the **vertical line**: corrected mass flow
is pinned to its design value `turb.mbD`.  The residual is:

    r = mb − turb.mbD

and its derivative with respect to `mb` is identically 1.  The design mass
flow `turb.mbD` is a constant anchor (no Newton-variable dependence), so its
contribution to the Newton Jacobian is zero.

This function is the canonical form used by the Newton driver to assemble
rows 2 and 3 of the 9×9 turbine residual system (HPT and LPT corrected-flow
constraints).  The full Jacobian column entries `a[i,j] = dr/d(newton_j)` are
obtained by chain-ruling `dr_dmb * dmb_d(newton_j) = 1 × mbXX_j` in the
Newton assembly.

## Arguments

| Argument | Unit  | Description                        |
|:---------|:------|:-----------------------------------|
| `turb`   | —     | `Turbine` component                |
| `mb`     | kg/s  | Off-design corrected mass flow     |

## Returns

`(r, dr_dmb)` where `r = mb - turb.mbD` and `dr_dmb = 1`.
"""
function turbine_mb_residual(turb::Turbine{T}, mb::T) where {T<:AbstractFloat}
    r      = mb - turb.mbD
    dr_dmb = one(T)
    return r, dr_dmb
end

# ---------------------------------------------------------------------------
# turbine_exit!  —  compute outlet FlowStation from inlet + work extraction
# ---------------------------------------------------------------------------

function turbine_exit!(
    st_out ::FlowStation{T},
    st_in  ::FlowStation{T},
    dh     ::T,
    ept    ::T,
) where {T<:AbstractFloat}
    epi = one(T) / ept

    # gas_delhd returns:
    #   pt, Tt, ht, s, cp, r, p_so, p_po, p_ep,
    #   p_ho, t_ho, h_ho, s_ho, p_dh,
    #   t_dh, h_dh, s_dh, p_al, t_al,
    #   h_al, s_al, cp_al, r_al
    res = gas_delhd(
        st_in.alpha, 5,
        st_in.pt, st_in.Tt, st_in.ht, st_in.st, st_in.cpt, st_in.Rt,
        dh, epi,
    )

    st_out.pt    = res[1]
    st_out.Tt    = res[2]
    st_out.ht    = res[3]
    st_out.st    = res[4]
    st_out.cpt   = res[5]
    st_out.Rt    = res[6]
    st_out.alpha = st_in.alpha   # composition unchanged through turbine

    return st_out
end
