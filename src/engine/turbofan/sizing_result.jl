"""
    SizingResult{T<:Real}

Named container for the outputs of `tfsize!`, parametric on the numeric type `T`.

Fields correspond 1-to-1 to the positional tuple elements previously returned by
`tfsize!`, in the same order.  Replacing the 50-variable positional destructure
with named field access eliminates silent mismatch between the return position and
the caller's variable name.

The type parameter `T` tracks the element type of all thermodynamic scalars and
the cooling arrays.  `T = Float64` in normal operation; `T = ForwardDiff.Dual{...}`
when the caller differentiates through `tfsize!` via forward-mode AD.  The
convergence flag `Lconv::Bool` is non-numeric and is excluded from `T`.

Station naming convention mirrors the `EngineState` flow-station names (st0, st2,
st25, st3, st4, st41, st45, st5, st8, st18, etc.).  The static-state group at
positions 89–94 (static temperature/velocity/pressure/cp/R/area at station 25) uses
the suffix `_25s` to distinguish them from the station-25 total-state group above.

See also: [`tfsize!`](@ref), [`tfcalc!`](@ref).
"""
struct SizingResult{T<:Real}
    # Cooling arrays
    epsrow ::Vector{T}   # cooling flow fractions (ncrowx = 4 rows)
    Tmrow  ::Vector{T}   # blade-row metal temperatures [K]

    # Thermodynamic performance scalars
    TSFC  ::T   # thrust-specific fuel consumption [kg/N/s]
    Fsp   ::T   # specific thrust [N/(kg/s)]
    hfuel ::T   # fuel specific enthalpy [J/kg]
    ff    ::T   # fuel-to-air ratio [—]
    mcore ::T   # core mass flow [kg/s]

    # Station 0 — freestream totals
    Tt0 ::T; ht0 ::T; pt0 ::T; cpt0 ::T; Rt0 ::T

    # Station 12 — fan-face outer totals
    Tt12 ::T; ht12 ::T; pt12 ::T; cpt12 ::T; Rt12 ::T

    # Station 2a — fan-face LPC totals
    Tt2a ::T; ht2a ::T; pt2a ::T; cpt2a ::T; Rt2a ::T

    # Station 2ac — pre-cooler exit totals
    Tt2ac ::T; ht2ac ::T; pt2ac ::T; cpt2ac ::T; Rt2ac ::T

    # Station 2 — fan-face totals
    Tt2 ::T; ht2 ::T; pt2 ::T; cpt2 ::T; Rt2 ::T

    # Station 13 — fan exit totals
    Tt13 ::T; ht13 ::T; pt13 ::T; cpt13 ::T; Rt13 ::T

    # Station 25 — LPC exit totals
    Tt2_5 ::T; ht2_5 ::T; pt2_5 ::T; cpt2_5 ::T; Rt2_5 ::T

    # Station 25c — inter-cooler exit totals
    Tt2_5c ::T; ht2_5c ::T; pt2_5c ::T; cpt2_5c ::T; Rt2_5c ::T

    # Station 3 — HPC exit totals
    Tt3 ::T; ht3 ::T; pt3 ::T; cpt3 ::T; Rt3 ::T

    # Station 4 — combustor exit totals (Tt4 is an INPUT; only the derived outputs returned)
    ht4 ::T; pt4 ::T; cpt4 ::T; Rt4 ::T

    # Station 41 — turbine inlet totals (after cooling mixing)
    Tt4_1 ::T; ht4_1 ::T; pt4_1 ::T; cpt4_1 ::T; Rt4_1 ::T

    # Station 45 — HPT exit totals
    Tt4_5 ::T; ht4_5 ::T; pt4_5 ::T; cpt4_5 ::T; Rt4_5 ::T

    # Station 5 — LPT exit totals
    Tt5 ::T; ht5 ::T; pt5 ::T; cpt5 ::T; Rt5 ::T

    # Station 8 — core nozzle throat totals
    Tt8 ::T; ht8 ::T; pt8 ::T; cpt8 ::T; Rt8 ::T

    # Station 18 — fan nozzle throat totals
    Tt18 ::T; ht18 ::T; pt18 ::T; cpt18 ::T; Rt18 ::T

    # Freestream velocity (updated for BLI)
    u0 ::T

    # Station 2 statics
    T2 ::T; u2 ::T; p2 ::T; cp2 ::T; R2 ::T; A2 ::T

    # Station 25 statics (tfsize! internal name: T2_5c/u2_5c/…; tfcalc! historically T2_5/u2_5/…)
    T2_5 ::T; u2_5 ::T; p2_5 ::T; cp2_5 ::T; R2_5 ::T; A2_5 ::T

    # Station 5/8 statics (core nozzle throat)
    T5 ::T; u5 ::T; p5 ::T; cp5 ::T; R5 ::T; A5 ::T

    # Station 6/9 statics (core nozzle exit)
    T6 ::T; u6 ::T; p6 ::T; cp6 ::T; R6 ::T; A6 ::T

    # Station 7/18 statics (fan nozzle throat)
    T7 ::T; u7 ::T; p7 ::T; cp7 ::T; R7 ::T; A7 ::T

    # Station 8/19 statics (fan nozzle exit)
    T8 ::T; u8 ::T; p8 ::T; cp8 ::T; R8 ::T; A8 ::T

    # Offtake discharge
    u9 ::T; A9 ::T

    # Component polytropic loss fractions
    epf ::T; eplc ::T; ephc ::T; epht ::T; eplt ::T

    # Component adiabatic efficiencies
    etaf ::T; etalc ::T; etahc ::T; etaht ::T; etalt ::T

    # Convergence flag
    Lconv ::Bool
end
