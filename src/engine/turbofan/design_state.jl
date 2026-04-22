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

## Cooling mixing constants
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
    pi_fan_des  ::T   # fan design total-pressure ratio
    pi_lpc_des ::T   # LPC design total-pressure ratio
    pi_hpc_des ::T   # HPC design total-pressure ratio
    pi_hpt_des ::T   # HPT design total-pressure ratio
    pi_lpt_des ::T   # LPT design total-pressure ratio

    # Map scalars — design-point corrected mass flows [kg/s · √K / kPa]
    mb_fan_des  ::T   # fan
    mb_lpc_des ::T   # LPC
    mb_hpc_des ::T   # HPC
    mb_hpt_des ::T   # HPT
    mb_lpt_des ::T   # LPT

    # Map scalars — design-point corrected spool speeds (normalised)
    Nb_fan_des  ::T   # fan / LP spool speed at design
    Nb_lpc_des ::T   # LPC design-point corrected speed
    Nb_hpc_des ::T   # HPC design-point corrected speed
    Nb_hpt_des ::T   # HPT design-point corrected speed
    Nb_lpt_des ::T   # LPT design-point corrected speed

    # -----------------------------------------------------------------------
    # Component flow areas [m²]
    # -----------------------------------------------------------------------
    A2    ::T   # fan-face area at station 2
    A25   ::T   # HPC-inlet area at station 25 (after intercooler)
    A8    ::T   # core nozzle throat area at station 8 (ARP755)
    A18   ::T   # fan nozzle throat area at station 18 (ARP755)

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

    # -----------------------------------------------------------------------
    # Component pressure ratios — frozen design inputs (iepid..iepitn)
    # -----------------------------------------------------------------------
    pid    ::T   # diffuser total-pressure ratio
    pib    ::T   # burner total-pressure ratio
    pifn   ::T   # fan nozzle total-pressure ratio
    pitn   ::T   # core nozzle total-pressure ratio

    # -----------------------------------------------------------------------
    # Polytropic efficiencies — frozen design inputs (ieepolf..ieepollt)
    # -----------------------------------------------------------------------
    epolf  ::T   # fan polytropic efficiency
    epollc ::T   # LPC polytropic efficiency
    epolhc ::T   # HPC polytropic efficiency
    epolht ::T   # HPT polytropic efficiency
    epollt ::T   # LPT polytropic efficiency

    # -----------------------------------------------------------------------
    # Fan map constants — frozen design inputs (iepifK, ieepfK)
    # -----------------------------------------------------------------------
    pifK   ::T   # fan pressure ratio at design corrected speed (FPR0)
    epfK   ::T   # fan polytropic efficiency factor (Kf)

    # -----------------------------------------------------------------------
    # Duct design Mach numbers — frozen inputs
    # -----------------------------------------------------------------------
    M2     ::T   # fan-face (station 2) design Mach number
    M25    ::T   # HPC-inlet (station 25) design Mach number

    # -----------------------------------------------------------------------
    # Spool mechanical losses — frozen inputs
    # -----------------------------------------------------------------------
    epsl   ::T   # LP spool mechanical loss fraction
    epsh   ::T   # HP spool mechanical loss fraction

    # -----------------------------------------------------------------------
    # Turbomachinery geometry — frozen TOML inputs (parg[igGearf/igHTRf/...])
    # -----------------------------------------------------------------------
    Gearf  ::T   # fan gear ratio (Nl/Nf); 1 = direct-drive
    HTRf   ::T   # fan hub-to-tip radius ratio
    HTRlc  ::T   # LPC hub-to-tip radius ratio
    HTRhc  ::T   # HPC hub-to-tip radius ratio

    # -----------------------------------------------------------------------
    # Combustion efficiency — frozen design input (ieetab)
    # -----------------------------------------------------------------------
    etab   ::T   # combustor adiabatic efficiency

    # -----------------------------------------------------------------------
    # Cooling design parameters — frozen inputs (iedTstrk..ietfilm)
    # -----------------------------------------------------------------------
    dTstrk  ::T   # hot-streak temperature allowance [K]
    Mtexit  ::T   # design Mach number at turbine blade exit
    StA     ::T   # Stanton number times blade area (St × A)
    efilm   ::T   # film cooling effectiveness
    tfilm   ::T   # film cooling thickness parameter

    # -----------------------------------------------------------------------
    # HPT cooling model constants — frozen design inputs (iefc0, iedehtdfc)
    # -----------------------------------------------------------------------
    fc0     ::T   # baseline cooling fraction (reference f_c for efficiency model)
    dehtdfc ::T   # HPT isentropic efficiency derivative w.r.t. cooling fraction
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
        z, z, z, z, z,   # pi_fan_des, pi_lpc_des, pi_hpc_des, pi_hpt_des, pi_lpt_des
        z, z, z, z, z,   # mb_fan_des, mb_lpc_des, mb_hpc_des, mb_hpt_des, mb_lpt_des
        z, z, z, z, z,   # Nb_fan_des, Nb_lpc_des, Nb_hpc_des, Nb_hpt_des, Nb_lpt_des
        z, z, z, z,      # A2, A25, A8, A18
        zv, zv,          # epsrow, Tmrow
        z,               # fc
        z, z,            # ruc, M4a
        z, z, z, z,      # pid, pib, pifn, pitn
        z, z, z, z, z,   # epolf, epollc, epolhc, epolht, epollt
        z, z,            # pifK, epfK
        z, z,            # M2, M25
        z, z,            # epsl, epsh
        z, z, z, z,     # Gearf, HTRf, HTRlc, HTRhc
        z,               # etab
        z, z, z, z, z,   # dTstrk, Mtexit, StA, efilm, tfilm
        z, z,            # fc0, dehtdfc
    )
end

"""
    DesignState()

Return a zero-initialised `DesignState{Float64}`.
"""
DesignState() = DesignState{Float64}()
