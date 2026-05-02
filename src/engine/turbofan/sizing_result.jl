"""
    SizingResult

Named container for the outputs of `tfsize!`.

Fields correspond 1-to-1 to the positional tuple elements previously returned by
`tfsize!`, in the same order.  Replacing the 50-variable positional destructure
with named field access eliminates silent mismatch between the return position and
the caller's variable name.

Station naming convention mirrors the `EngineState` flow-station names (st0, st2,
st25, st3, st4, st41, st45, st5, st8, st18, etc.).  The static-state group at
positions 89–94 (static temperature/velocity/pressure/cp/R/area at station 25) uses
the suffix `_25s` to distinguish them from the station-25 total-state group above.

See also: [`tfsize!`](@ref), [`tfcalc!`](@ref).
"""
struct SizingResult
    # Cooling arrays
    epsrow ::Vector{Float64}   # cooling flow fractions (ncrowx = 4 rows)
    Tmrow  ::Vector{Float64}   # blade-row metal temperatures [K]

    # Thermodynamic performance scalars
    TSFC  ::Float64   # thrust-specific fuel consumption [kg/N/s]
    Fsp   ::Float64   # specific thrust [N/(kg/s)]
    hfuel ::Float64   # fuel specific enthalpy [J/kg]
    ff    ::Float64   # fuel-to-air ratio [—]
    mcore ::Float64   # core mass flow [kg/s]

    # Station 0 — freestream totals
    Tt0 ::Float64; ht0 ::Float64; pt0 ::Float64; cpt0 ::Float64; Rt0 ::Float64

    # Station 12 — fan-face outer totals
    Tt12 ::Float64; ht12 ::Float64; pt12 ::Float64; cpt12 ::Float64; Rt12 ::Float64

    # Station 2a — fan-face LPC totals
    Tt2a ::Float64; ht2a ::Float64; pt2a ::Float64; cpt2a ::Float64; Rt2a ::Float64

    # Station 2ac — pre-cooler exit totals
    Tt2ac ::Float64; ht2ac ::Float64; pt2ac ::Float64; cpt2ac ::Float64; Rt2ac ::Float64

    # Station 2 — fan-face totals
    Tt2 ::Float64; ht2 ::Float64; pt2 ::Float64; cpt2 ::Float64; Rt2 ::Float64

    # Station 13 — fan exit totals
    Tt13 ::Float64; ht13 ::Float64; pt13 ::Float64; cpt13 ::Float64; Rt13 ::Float64

    # Station 25 — LPC exit totals
    Tt2_5 ::Float64; ht2_5 ::Float64; pt2_5 ::Float64; cpt2_5 ::Float64; Rt2_5 ::Float64

    # Station 25c — inter-cooler exit totals
    Tt2_5c ::Float64; ht2_5c ::Float64; pt2_5c ::Float64; cpt2_5c ::Float64; Rt2_5c ::Float64

    # Station 3 — HPC exit totals
    Tt3 ::Float64; ht3 ::Float64; pt3 ::Float64; cpt3 ::Float64; Rt3 ::Float64

    # Station 4 — combustor exit totals (Tt4 is an INPUT; only the derived outputs returned)
    ht4 ::Float64; pt4 ::Float64; cpt4 ::Float64; Rt4 ::Float64

    # Station 41 — turbine inlet totals (after cooling mixing)
    Tt4_1 ::Float64; ht4_1 ::Float64; pt4_1 ::Float64; cpt4_1 ::Float64; Rt4_1 ::Float64

    # Station 45 — HPT exit totals
    Tt4_5 ::Float64; ht4_5 ::Float64; pt4_5 ::Float64; cpt4_5 ::Float64; Rt4_5 ::Float64

    # Station 5 — LPT exit totals
    Tt5 ::Float64; ht5 ::Float64; pt5 ::Float64; cpt5 ::Float64; Rt5 ::Float64

    # Station 8 — core nozzle throat totals
    Tt8 ::Float64; ht8 ::Float64; pt8 ::Float64; cpt8 ::Float64; Rt8 ::Float64

    # Station 18 — fan nozzle throat totals
    Tt18 ::Float64; ht18 ::Float64; pt18 ::Float64; cpt18 ::Float64; Rt18 ::Float64

    # Freestream velocity (updated for BLI)
    u0 ::Float64

    # Station 2 statics
    T2 ::Float64; u2 ::Float64; p2 ::Float64; cp2 ::Float64; R2 ::Float64; A2 ::Float64

    # Station 25 statics (tfsize! internal name: T2_5c/u2_5c/…; tfcalc! historically T2_5/u2_5/…)
    T2_5 ::Float64; u2_5 ::Float64; p2_5 ::Float64; cp2_5 ::Float64; R2_5 ::Float64; A2_5 ::Float64

    # Station 5/8 statics (core nozzle throat)
    T5 ::Float64; u5 ::Float64; p5 ::Float64; cp5 ::Float64; R5 ::Float64; A5 ::Float64

    # Station 6/9 statics (core nozzle exit)
    T6 ::Float64; u6 ::Float64; p6 ::Float64; cp6 ::Float64; R6 ::Float64; A6 ::Float64

    # Station 7/18 statics (fan nozzle throat)
    T7 ::Float64; u7 ::Float64; p7 ::Float64; cp7 ::Float64; R7 ::Float64; A7 ::Float64

    # Station 8/19 statics (fan nozzle exit)
    T8 ::Float64; u8 ::Float64; p8 ::Float64; cp8 ::Float64; R8 ::Float64; A8 ::Float64

    # Offtake discharge
    u9 ::Float64; A9 ::Float64

    # Component polytropic loss fractions
    epf ::Float64; eplc ::Float64; ephc ::Float64; epht ::Float64; eplt ::Float64

    # Component adiabatic efficiencies
    etaf ::Float64; etalc ::Float64; etahc ::Float64; etaht ::Float64; etalt ::Float64

    # Convergence flag
    Lconv ::Bool
end
