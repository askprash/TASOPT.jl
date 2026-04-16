"""
Inlet component type for a TASOPT turbofan.

Models the intake duct from the far-field (station 0) to the fan/core inlet
face (stations 2 and 19), including optional boundary-layer-ingestion (BLI)
entropy mixing.
"""

# ---------------------------------------------------------------------------
# Inlet struct
# ---------------------------------------------------------------------------

"""
    Inlet{T<:AbstractFloat}

Parametric inlet component holding pressure-recovery and BLI (boundary layer
ingestion) parameters.

Parametric in the numeric type `T` so that forward-mode AD and other dual-
number types flow through without requiring specialised containers.

## Fields

| Field               | Unit | Description                                            |
|:--------------------|:-----|:-------------------------------------------------------|
| `pid`               | —    | Diffuser total-pressure ratio `pt18 / pt0`, typically 0.97–1.0 |
| `Kinl`              | W    | BLI kinetic-energy defect ingested from the airframe boundary layer; zero for clean (non-BLI) inlets |
| `Phiinl`            | W    | BLI ingested dissipation (mechanical power); used for the airframe mechanical-energy balance; zero for clean inlets |
| `eng_has_BLI_cores` | —    | If `true`, both fan and core streams see the BLI entropy defect; if `false`, only the fan stream is affected |

## Physics summary

The inlet is split into two stages:

1. **Adiabatic diffuser (0 → 18):** Total enthalpy and temperature are
   unchanged; total pressure drops by `pid`:
   ```
   pt18 = pt0 × pid,   Tt18 = Tt0
   ```

2. **BLI entropy mixing (18 → 2/19):** The airframe boundary layer, whose
   kinetic-energy defect is `Kinl`, mixes with the incoming stream and
   causes an additional entropy rise that further reduces total pressure:
   ```
   sbfan  = Kinl × γ₀ / (ṁ_mix × a₂²)
   pt2    = pt18 × exp(−sbfan)
   pt19   = pt18 × exp(−sbcore)    # sbcore = sbfan if eng_has_BLI_cores, else 0
   Tt2    = Tt19 = Tt18            # adiabatic BLI
   ```
   where `ṁ_mix` is the BLI-ingested mass flow (approximated from corrected
   flows) and `a₂` is the speed of sound at the fan-face Mach number.
   When `Kinl = 0` or `u0 = 0` (static condition), the BLI mixing step is
   a no-op and stations 2/19 equal station 18.
"""
mutable struct Inlet{T<:AbstractFloat}
    pid              ::T    # diffuser total-pressure ratio  pt18/pt0   [—]
    Kinl             ::T    # BLI kinetic-energy defect                 [W]
    Phiinl           ::T    # BLI ingested dissipation                  [W]
    eng_has_BLI_cores::Bool # core stream also sees BLI defect?
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    Inlet(pid::T; Kinl=zero(T), Phiinl=zero(T), eng_has_BLI_cores=false)

Convenience constructor for a clean (BLI-free) inlet with pressure-recovery
ratio `pid`.  BLI parameters default to zero / false.
"""
function Inlet(
    pid              ::T;
    Kinl             ::T    = zero(T),
    Phiinl           ::T    = zero(T),
    eng_has_BLI_cores::Bool = false,
) where {T<:AbstractFloat}
    Inlet{T}(pid, Kinl, Phiinl, eng_has_BLI_cores)
end

# ---------------------------------------------------------------------------
# inlet_diffuser! — adiabatic diffuser, station 0 → 18
# ---------------------------------------------------------------------------

"""
    inlet_diffuser!(st18, st0, inlet) -> st18

Compute station 18 (inlet face after the diffuser) from the freestream station
`st0` and the inlet's pressure-recovery ratio `inlet.pid`.

**Adiabatic duct assumption:** total enthalpy, temperature, entropy complement,
specific heats, gas constant, and species composition are all copied unchanged
from `st0`.  Only total pressure changes:

    pt18 = pt0 × inlet.pid

Static-state fields (`Ts`, `ps`, `hs`, `u`, etc.) are **not** set here because
they depend on the local Mach number, which is determined by downstream flow
matching; callers must compute the static state separately via
`set_static_from_M!` once the Mach number is known.

Returns `st18` (mutated in place).
"""
function inlet_diffuser!(
    st18 ::FlowStation{T},
    st0  ::FlowStation{T},
    inlet::Inlet,
) where {T<:AbstractFloat}
    # Adiabatic: total enthalpy/temperature/entropy unchanged
    st18.Tt    = st0.Tt
    st18.ht    = st0.ht
    st18.cpt   = st0.cpt
    st18.Rt    = st0.Rt
    st18.st    = st0.st
    st18.alpha = st0.alpha

    # Pressure recovery
    st18.pt = st0.pt * inlet.pid

    return st18
end

# ---------------------------------------------------------------------------
# inlet_bli_mixing! — BLI entropy mixing, station 18 → 2 / 19
# ---------------------------------------------------------------------------

"""
    inlet_bli_mixing!(st2, st19, st18, st0, inlet,
                      mf, ml, M2, at0, gam0, Tref, pref)
      -> (; sbfan, sbcore)

Apply BLI (boundary-layer ingestion) entropy mixing from station 18 to
stations 2 (fan face) and 19 (core inlet).

## Arguments

| Argument  | Unit   | Description                                                  |
|:----------|:-------|:-------------------------------------------------------------|
| `st2`     | —      | Fan-face station (mutated)                                   |
| `st19`    | —      | LPC-inlet station (mutated)                                  |
| `st18`    | —      | Inlet-face station (read-only)                               |
| `st0`     | —      | Freestream station; used for freestream `pt0`, `Tt0` in BLI reference calculation |
| `inlet`   | —      | `Inlet` component holding `Kinl` and `eng_has_BLI_cores`    |
| `mf`      | —      | Dimensionless corrected fan mass flow (Newton iterate, ≈1 at design) |
| `ml`      | —      | Dimensionless corrected LPC mass flow (Newton iterate)       |
| `M2`      | —      | Fan-inlet Mach number (Newton iterate)                       |
| `at0`     | m/s    | Stagnation speed of sound at freestream                      |
| `gam0`    | —      | Ratio of specific heats at freestream                        |
| `Tref`    | K      | Reference temperature for corrected mass-flow normalisation  |
| `pref`    | Pa     | Reference pressure for corrected mass-flow normalisation     |

## Physics

The BLI entropy boost is modelled as a mass-averaged entropy increase over the
ingested flow:

    a₂²    = at0² / (1 + 0.5·(γ₀−1)·M₂²)
    ṁ_mix  = mf·√(Tref/Tt0)·(pt0/pref)          [fan-only case]
           + ml·√(Tref/Tt0)·(pt0/pref)          [+ core term if BLI_cores]
    sbfan  = Kinl·γ₀ / (ṁ_mix·a₂²)
    pt2    = pt18·exp(−sbfan)
    pt19   = pt18·exp(−sbcore)    # sbcore = sbfan or 0

where the approximation `pt0 ≈ pt2`, `Tt0 ≈ Tt2` (used in `ṁ_mix`) avoids the
circular dependency between pt2 and sbfan (same approximation used in
`tfoper!`).

When `Kinl = 0` or `at0 = 0` (static condition), all BLI terms vanish and
stations 2/19 receive the state of station 18 directly.

Total temperature is conserved (adiabatic BLI assumption):
`Tt2 = Tt19 = Tt18`.

## Returns

A named tuple `(; sbfan, sbcore)` with the entropy-boost factors for both
streams.  These are useful for analytic Jacobian assembly in the Newton driver
(tasopt-j9l.31).
"""
function inlet_bli_mixing!(
    st2  ::FlowStation{T},
    st19 ::FlowStation{T},
    st18 ::FlowStation{T},
    st0  ::FlowStation{T},
    inlet::Inlet,
    mf   ::T, ml::T, M2::T,
    at0  ::T, gam0::T, Tref::T, pref::T,
) where {T<:AbstractFloat}
    # Adiabatic BLI: total enthalpy/temperature/entropy unchanged
    for st_out in (st2, st19)
        st_out.Tt    = st18.Tt
        st_out.ht    = st18.ht
        st_out.cpt   = st18.cpt
        st_out.Rt    = st18.Rt
        st_out.st    = st18.st
        st_out.alpha = st18.alpha
    end

    z = zero(T)

    # Static condition (u0=0) or no BLI: entropy defect is zero
    if inlet.Kinl == z || at0 == z
        st2.pt  = st18.pt
        st19.pt = st18.pt
        return (; sbfan=z, sbcore=z)
    end

    # Speed of sound squared at fan-face Mach number
    a2sq = at0^2 / (one(T) + T(0.5) * (gam0 - one(T)) * M2^2)

    # Freestream reference factors (avoiding circular dependency on pt2/pt19)
    pt0  = st0.pt
    Tt0  = st0.Tt
    ref_scale = sqrt(Tref / Tt0) * (pt0 / pref)

    local sbfan, sbcore
    if inlet.eng_has_BLI_cores
        # BL mixes with fan + core flow
        mmix   = (mf + ml) * ref_scale
        sbfan  = inlet.Kinl * gam0 / (mmix * a2sq)
        sbcore = sbfan
    else
        # BL mixes with fan flow only
        mmix   = mf * ref_scale
        sbfan  = inlet.Kinl * gam0 / (mmix * a2sq)
        sbcore = z
    end

    st2.pt  = st18.pt * exp(-sbfan)
    st19.pt = st18.pt * exp(-sbcore)

    return (; sbfan, sbcore)
end
