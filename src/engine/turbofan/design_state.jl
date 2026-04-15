"""
Frozen design-point state container for a TASOPT turbofan.
"""

using StaticArrays

# ---------------------------------------------------------------------------
# DesignState — frozen scalars set during tfsize!, read during tfoper!
# ---------------------------------------------------------------------------

"""
    DesignState{T<:AbstractFloat}

Immutable (frozen) design-point scalars for a TASOPT turbofan engine.
These quantities are computed once during on-design sizing (`tfsize!`) and
held fixed throughout all subsequent off-design evaluations (`tfoper!`).

Parametric in the numeric type `T` so that forward-mode AD (ForwardDiff)
and other dual-number types flow through without requiring specialised
containers.

## Map scalars
Component performance maps are normalised relative to the design-point
corrected speeds (`Nb*D`), corrected mass flows (`mb*D`), and pressure
ratios (`pi*D`).  Subscript conventions: f = fan, lc = LPC, hc = HPC,
ht = HPT, lt = LPT.

## Component areas
Flow areas at key cross-sections, sized to pass the design corrected flow at
the design Mach number.

## Cooling
A maximum of `ncrowx = 4` cooled blade rows are supported (this is a fixed
dimension of the current cooling model; no run-time variant exists).
Depending on `CoolingOpt`, either the bypass ratios (`epsrow`) or the metal
temperatures (`Tmrow`) are specified inputs and the other is computed during
sizing; both are stored here after the design is frozen.

## OQ-4 decision (aerospace + programmer, 2026-04-14)
`ruc` (cooling-flow velocity ratio at the mixing plane) and `M4a`
(prescribed Mach number at the start-of-mixing / cooling-flow outlet,
station 4a) are frozen design constants — not thermodynamic station
variables — with zero Newton partial derivatives.  They live here, not in a
`FlowStation`.
"""
mutable struct DesignState{T<:AbstractFloat}
    # -----------------------------------------------------------------------
    # Map scalars — design-point pressure ratios (dimensionless)
    # -----------------------------------------------------------------------
    pifD  ::T   # fan design total-pressure ratio
    pilcD ::T   # LPC design total-pressure ratio
    pihcD ::T   # HPC design total-pressure ratio
    pihtD ::T   # HPT design total-pressure ratio
    piltD ::T   # LPT design total-pressure ratio

    # Map scalars — design-point corrected mass flows [kg/s · √K / kPa]
    mbfD  ::T   # fan
    mblcD ::T   # LPC
    mbhcD ::T   # HPC
    mbhtD ::T   # HPT
    mbltD ::T   # LPT

    # Map scalars — design-point corrected spool speeds (normalised)
    NbfD  ::T   # fan / LP spool speed at design
    NblcD ::T   # LPC design-point corrected speed
    NbhcD ::T   # HPC design-point corrected speed
    NbhtD ::T   # HPT design-point corrected speed
    NbltD ::T   # LPT design-point corrected speed

    # -----------------------------------------------------------------------
    # Component flow areas [m²]
    # -----------------------------------------------------------------------
    A2    ::T   # fan-face area at station 2
    A25   ::T   # HPC-inlet area at station 25 (after intercooler)
    A5    ::T   # core nozzle throat area at station 5
    A7    ::T   # fan nozzle throat area at station 7

    # -----------------------------------------------------------------------
    # Cooling (ncrowx = 4 blade rows — fixed compile-time dimension)
    # -----------------------------------------------------------------------
    # Cooling-flow bypass ratio for each blade row:
    #   epsrow[i] = ṁ_cool,i / ṁ_core,design
    epsrow ::SVector{4,T}
    # Blade metal temperature for each blade row [K]
    Tmrow  ::SVector{4,T}
    # Total cooling mass-flow fraction (normalised by design core flow)
    fc     ::T

    # -----------------------------------------------------------------------
    # Cooling mixing constants (OQ-4 resolved: frozen design inputs)
    # -----------------------------------------------------------------------
    # Cooling-flow velocity ratio u_coolant / u_mainstream at station 4a
    ruc    ::T
    # Prescribed Mach number at the start-of-mixing station 4a
    # (bridges the momentum-weighted mixing from station 4 → 41)
    M4a    ::T
end

"""
    DesignState{T}() where {T<:AbstractFloat}

Return a zero-initialised `DesignState` with numeric type `T`.
All scalar fields are set to `zero(T)`; vector fields to zero `SVector`s.
"""
function DesignState{T}() where {T<:AbstractFloat}
    z = zero(T)
    zv = @SVector zeros(T, 4)
    DesignState{T}(
        z, z, z, z, z,   # pifD, pilcD, pihcD, pihtD, piltD
        z, z, z, z, z,   # mbfD, mblcD, mbhcD, mbhtD, mbltD
        z, z, z, z, z,   # NbfD, NblcD, NbhcD, NbhtD, NbltD
        z, z, z, z,      # A2, A25, A5, A7
        zv, zv,          # epsrow, Tmrow
        z,               # fc
        z, z,            # ruc, M4a
    )
end

"""
    DesignState()

Return a zero-initialised `DesignState{Float64}`.
"""
DesignState() = DesignState{Float64}()
