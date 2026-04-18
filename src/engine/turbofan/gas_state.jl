"""
Gas-state container for a single thermodynamic control volume.
"""

using StaticArrays

# ---------------------------------------------------------------------------
# GasState — total + static thermo state and 5-species mass fractions
# ---------------------------------------------------------------------------

"""
    GasState{T<:AbstractFloat}

Mutable container for the complete thermodynamic state at a single flow
station in the TASOPT thermally-perfect gas model.

Parametric in the numeric type `T` so that forward-mode AD (ForwardDiff)
and other dual-number types flow through without requiring specialised
containers.

## Gas model

Five-species thermally-perfect mixture: N₂ (i=1), O₂ (i=2), CO₂ (i=3),
H₂O (i=4), Ar (i=5).  Species composition is stored as an `SVector{5,T}`
of mass fractions.  Pre-combustion air uses the conventional `alpha` vector;
post-combustion products use `lambda`.  Both representations have exactly
five components and live in the same field — callers are responsible for
tracking which interpretation applies.

## Notes

The legacy `gascalc.jl` interface passes `n` (number of species) as a
runtime `Int` argument at every call site (≥50 call sites across 13
functions).  Changing all those signatures in one step spans Phase 2/3/4
files and creates unacceptable merge risk at this stage.

Resolution — option (b): **keep `n` as `Int` in the legacy interface** and
add `@assert n == 5` in debug builds at each call site, so any inadvertent
mismatch is caught immediately.

## Fields

### Total (stagnation) state
| Field  | Unit        | Description                                         |
|:-------|:------------|:----------------------------------------------------|
| `Tt`   | K           | Total temperature                                   |
| `ht`   | J/kg        | Total complete enthalpy (includes heat of formation)|
| `pt`   | Pa          | Total pressure                                      |
| `cpt`  | J/(kg·K)    | Specific heat at stagnation temperature (= dh/dT)   |
| `Rt`   | J/(kg·K)    | Gas constant at stagnation conditions               |
| `st`   | J/(kg·K)    | Entropy-complement function s[Tt] at stagnation     |

### Static state
| Field  | Unit        | Description                                         |
|:-------|:------------|:----------------------------------------------------|
| `Ts`   | K           | Static temperature                                  |
| `ps`   | Pa          | Static pressure                                     |
| `hs`   | J/kg        | Static enthalpy                                     |
| `ss`   | J/(kg·K)    | Entropy-complement function s[T]                    |
| `cps`  | J/(kg·K)    | Specific heat at static temperature                 |
| `Rs`   | J/(kg·K)    | Gas constant at static conditions                   |
| `u`    | m/s         | Flow velocity                                       |

### Species composition
| Field   | —           | Description                                         |
|:--------|:------------|:----------------------------------------------------|
| `alpha` | kg/kg       | Mass fractions for 5 species (N₂, O₂, CO₂, H₂O, Ar)|
"""
mutable struct GasState{T<:AbstractFloat}
    # -----------------------------------------------------------------------
    # Total (stagnation) thermodynamic state
    # -----------------------------------------------------------------------
    Tt  ::T   # total temperature                            [K]
    ht  ::T   # total complete enthalpy (incl. heat of form) [J/kg]
    pt  ::T   # total pressure                               [Pa]
    cpt ::T   # specific heat at stagnation temperature      [J/(kg·K)]
    Rt  ::T   # gas constant at stagnation conditions        [J/(kg·K)]
    st  ::T   # entropy-complement s[Tt] at stagnation       [J/(kg·K)]

    # -----------------------------------------------------------------------
    # Static thermodynamic state
    # -----------------------------------------------------------------------
    Ts  ::T   # static temperature                           [K]
    ps  ::T   # static pressure                              [Pa]
    hs  ::T   # static enthalpy                              [J/kg]
    ss  ::T   # entropy-complement function s[T]             [J/(kg·K)]
    cps ::T   # specific heat at static temperature          [J/(kg·K)]
    Rs  ::T   # gas constant at static conditions            [J/(kg·K)]
    u   ::T   # flow velocity                                [m/s]

    # -----------------------------------------------------------------------
    # Species composition — 5 species, fixed at compile time
    #
    # OQ-5 (see docstring): n is NOT a field; n==5 is a compile-time
    # invariant enforced by the SVector length.  Legacy callers that pass n
    # as a runtime Int should @assert n == 5 when in debug mode.
    # -----------------------------------------------------------------------
    alpha ::SVector{5,T}   # mass fractions: [N₂, O₂, CO₂, H₂O, Ar]
end

"""
    GasState{T}() where {T<:AbstractFloat}

Return a zero-initialised `GasState` with numeric type `T`.
All scalar fields are `zero(T)`; `alpha` is a zero `SVector{5,T}`.
"""
function GasState{T}() where {T<:AbstractFloat}
    z  = zero(T)
    za = @SVector zeros(T, 5)
    GasState{T}(
        z, z, z, z, z, z,     # Tt, ht, pt, cpt, Rt, st
        z, z, z, z, z, z, z,  # Ts, ps, hs, ss, cps, Rs, u
        za,                    # alpha
    )
end

"""
    GasState()

Return a zero-initialised `GasState{Float64}`.
"""
GasState() = GasState{Float64}()

"""
    GasState{T}(Tt, ht, pt, cpt, Rt, alpha) where {T<:AbstractFloat}

Construct a `GasState` with explicit total-state fields and species
composition; all static-state fields default to `zero(T)`.

Useful when only stagnation conditions are known (e.g. directly after a
combustor model) and static quantities will be filled in later.
"""
function GasState{T}(Tt::T, ht::T, pt::T, cpt::T, Rt::T,
                     alpha::SVector{5,T}) where {T<:AbstractFloat}
    z = zero(T)
    GasState{T}(
        Tt, ht, pt, cpt, Rt,  # total state
        z,                     # st (entropy-complement at stagnation, not yet computed)
        z, z, z, z, z, z, z,  # static state (zeroed)
        alpha,
    )
end
