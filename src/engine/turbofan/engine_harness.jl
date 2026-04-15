"""
engine_harness.jl — Engine-standalone single-point runner.

Provides `run_engine_design_point`, which accepts an aircraft object (or a
path to an input TOML file), initialises the design-point ambient conditions,
and runs the turbofan design-point sizing routine (`tfsize!` via `tfwrap!`)
for one mission point — **without** iterating the full aircraft
structural/aerodynamic weight-convergence loop.

Also provides `pare_to_engine_state!`, the public mapping function that reads
a single `pare` column into a typed `EngineState`.
"""

# ---------------------------------------------------------------------------
# Internal helpers — station-level field population from a pare slice
# ---------------------------------------------------------------------------

"""
    _fill_total!(fs, pare, iTt, iht, ipt, icpt, iRt)

Write the five total-state scalars from a `pare` slice into `fs`.
"""
@inline function _fill_total!(fs::FlowStation, pare,
                               iTt::Int, iht::Int, ipt::Int, icpt::Int, iRt::Int)
    fs.Tt  = pare[iTt]
    fs.ht  = pare[iht]
    fs.pt  = pare[ipt]
    fs.cpt = pare[icpt]
    fs.Rt  = pare[iRt]
end

"""
    _fill_static!(fs, pare, ip, iT, iR, icp, iu)

Write the five static-state scalars (p, T, R, cp, u) from a `pare` slice
into `fs`.  Static enthalpy `hs` and entropy complement `ss` are not in
`pare`; they remain at their prior value.
"""
@inline function _fill_static!(fs::FlowStation, pare,
                                ip::Int, iT::Int, iR::Int, icp::Int, iu::Int)
    fs.ps  = pare[ip]
    fs.Ts  = pare[iT]
    fs.Rs  = pare[iR]
    fs.cps = pare[icp]
    fs.u   = pare[iu]
end

# ---------------------------------------------------------------------------
# pare_to_engine_state!
# ---------------------------------------------------------------------------

"""
    pare_to_engine_state!(eng, pare) -> eng

Populate a (typically zero-initialised) `EngineState` from a single column
of the `pare` array — `view(ac.pare, :, ip, imission)` — after a successful
`tfcalc!` / `tfsize!` call.

Only fields that the legacy `pare` array actually stores are written:

- **Total state** (Tt, ht, pt, cpt, Rt): stations 0, 18, 19, 2, 21, 25, 3,
  4, 41, 45, 49, 5, 7; partial for station 9 (Tt, pt only).
- **Static state** (ps, Ts, Rs, cps, u): stations 2, 25, 5, 6, 7, 8.
  Stations 6 and 8 have static state but no total state in `pare`.
- **Area** (A): stations 2, 25, 5, 6, 7, 8, 9.
- **Mass flow** (mdot): station 2 (`mcore`, the core mass flow).
- **Ambient scalars**: `M0`, `T0`, `p0`, `a0`.

Fields **not** in `pare` — total state for stations 19c, 25c, 4a, 49c; static
enthalpy `hs` and entropy complement `ss` at all stations — remain at their
prior value (zero for a fresh `EngineState()`).
"""
function pare_to_engine_state!(eng::EngineState, pare)
    # -----------------------------------------------------------------------
    # Ambient scalars
    # -----------------------------------------------------------------------
    eng.M0 = pare[ieM0]
    eng.T0 = pare[ieT0]
    eng.p0 = pare[iep0]
    eng.a0 = pare[iea0]

    # -----------------------------------------------------------------------
    # Station 0 — freestream
    # -----------------------------------------------------------------------
    _fill_total!(eng.st0, pare, ieTt0, ieht0, iept0, iecpt0, ieRt0)
    eng.st0.u = pare[ieu0]

    # -----------------------------------------------------------------------
    # Station 18 — FanFaceOuter
    # -----------------------------------------------------------------------
    _fill_total!(eng.st18, pare, ieTt18, ieht18, iept18, iecpt18, ieRt18)

    # -----------------------------------------------------------------------
    # Station 19 — FanFaceLPC
    # -----------------------------------------------------------------------
    _fill_total!(eng.st19, pare, ieTt19, ieht19, iept19, iecpt19, ieRt19)

    # Station 19c (PreCoolerOut) — NOT in pare; remains zero.

    # -----------------------------------------------------------------------
    # Station 2 — FanFaceFan  (total + static + area + mdot)
    # -----------------------------------------------------------------------
    _fill_total!(eng.st2, pare, ieTt2, ieht2, iept2, iecpt2, ieRt2)
    _fill_static!(eng.st2, pare, iep2, ieT2, ieR2, iecp2, ieu2)
    eng.st2.A    = pare[ieA2]
    eng.st2.mdot = pare[iemcore]

    # -----------------------------------------------------------------------
    # Station 21 — FanExit
    # -----------------------------------------------------------------------
    _fill_total!(eng.st21, pare, ieTt21, ieht21, iept21, iecpt21, ieRt21)

    # -----------------------------------------------------------------------
    # Station 25 — LPCExit  (total + static + area)
    # -----------------------------------------------------------------------
    _fill_total!(eng.st25, pare, ieTt25, ieht25, iept25, iecpt25, ieRt25)
    _fill_static!(eng.st25, pare, iep25, ieT25, ieR25, iecp25, ieu25)
    eng.st25.A = pare[ieA25]

    # Station 25c (InterCoolerOut) — NOT in pare; remains zero.

    # -----------------------------------------------------------------------
    # Station 3 — HPCExit
    # -----------------------------------------------------------------------
    _fill_total!(eng.st3, pare, ieTt3, ieht3, iept3, iecpt3, ieRt3)

    # -----------------------------------------------------------------------
    # Station 4 — CombustorExit
    # Tt4 is an INPUT to tfsize! (not written back by tfcalc!);
    # ht4 / pt4 / cpt4 / Rt4 are computed outputs.
    # -----------------------------------------------------------------------
    eng.st4.Tt  = pare[ieTt4]
    eng.st4.ht  = pare[ieht4]
    eng.st4.pt  = pare[iept4]
    eng.st4.cpt = pare[iecpt4]
    eng.st4.Rt  = pare[ieRt4]

    # Station 4a (CoolMixInlet) — NOT in pare; remains zero.

    # -----------------------------------------------------------------------
    # Station 41 — TurbineInlet
    # -----------------------------------------------------------------------
    _fill_total!(eng.st41, pare, ieTt41, ieht41, iept41, iecpt41, ieRt41)

    # -----------------------------------------------------------------------
    # Station 45 — HPTExit
    # -----------------------------------------------------------------------
    _fill_total!(eng.st45, pare, ieTt45, ieht45, iept45, iecpt45, ieRt45)

    # -----------------------------------------------------------------------
    # Station 49 — LPTExit
    # -----------------------------------------------------------------------
    _fill_total!(eng.st49, pare, ieTt49, ieht49, iept49, iecpt49, ieRt49)

    # Station 49c (RegenCoolerOut) — NOT in pare; remains zero.

    # -----------------------------------------------------------------------
    # Station 5 — CoreNozzle  (total + static + area)
    # -----------------------------------------------------------------------
    _fill_total!(eng.st5, pare, ieTt5, ieht5, iept5, iecpt5, ieRt5)
    _fill_static!(eng.st5, pare, iep5, ieT5, ieR5, iecp5, ieu5)
    eng.st5.A = pare[ieA5]

    # -----------------------------------------------------------------------
    # Station 6 — CoreNozzleExit  (static ONLY — no total in pare)
    # -----------------------------------------------------------------------
    _fill_static!(eng.st6, pare, iep6, ieT6, ieR6, iecp6, ieu6)
    eng.st6.A = pare[ieA6]

    # -----------------------------------------------------------------------
    # Station 7 — FanNozzle  (total + static + area)
    # -----------------------------------------------------------------------
    _fill_total!(eng.st7, pare, ieTt7, ieht7, iept7, iecpt7, ieRt7)
    _fill_static!(eng.st7, pare, iep7, ieT7, ieR7, iecp7, ieu7)
    eng.st7.A = pare[ieA7]

    # -----------------------------------------------------------------------
    # Station 8 — FanNozzleExit  (static ONLY — no total in pare)
    # -----------------------------------------------------------------------
    _fill_static!(eng.st8, pare, iep8, ieT8, ieR8, iecp8, ieu8)
    eng.st8.A = pare[ieA8]

    # -----------------------------------------------------------------------
    # Station 9 — OfftakeDisch  (Tt, pt, u, A; no ht/cpt/Rt in pare)
    # -----------------------------------------------------------------------
    eng.st9.Tt = pare[ieTt9]
    eng.st9.pt = pare[iept9]
    eng.st9.u  = pare[ieu9]
    eng.st9.A  = pare[ieA9]

    return eng
end

# ---------------------------------------------------------------------------
# run_engine_design_point
# ---------------------------------------------------------------------------

"""
    run_engine_design_point(ac; imission=1, ip=ipcruise1) -> EngineState

Run the turbofan design-point thermodynamic sizing for one mission point and
return the result as a populated `EngineState`.

The full structural/aerodynamic aircraft weight-convergence loop is **not**
run — only the engine thermodynamics (`tfsize!` via `tfcalc!`/`tfwrap!`).
Ambient conditions are computed from the ISA standard atmosphere at the
design altitude and Mach number declared in the input file.

## Arguments
- `ac`: an `aircraft` object loaded with `load_default_model()` or
  `read_aircraft_model()`.  `pare[ieFe, ip, imission]` must hold a meaningful
  design thrust before calling this function; the aircraft should be sized
  with `size_aircraft!` first, or `pare[ieFe]` set manually.
- `imission::Int`: mission index (default `1`).
- `ip::Int`: mission-point index (default `ipcruise1`).

## Returns
An `EngineState{Float64}` populated from the converged design-point `pare`
column at `ip`.  Station fields not stored in `pare` (total state for 19c,
25c, 4a, 49c; `hs` and `ss` for all stations) have value zero.

## Notes
For a fully-sized aircraft, `run_engine_design_point` reproduces the same
thermodynamic result as the final iteration of `size_aircraft!`.

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
size_aircraft!(ac; printiter=false)   # populates pare[ieFe]
eng = TASOPT.engine.run_engine_design_point(ac)
TASOPT.engine.dump_stations(eng)
```
"""
function run_engine_design_point(ac; imission::Int=1, ip::Int=ipcruise1)

    # -----------------------------------------------------------------------
    # Initialise atmosphere ISA offset — mirrors step 1 of _size_aircraft!
    # so that the standard-atmosphere call uses the correct ΔT.
    # -----------------------------------------------------------------------
    altTO = ac.parm[imaltTO]
    T_std = atmos(altTO).T
    ac.parm[imDeltaTatm] = ac.parm[imT0TO] - T_std

    # -----------------------------------------------------------------------
    # Set T0, p0, a0, M0, u0, ρ0, μ0 in pare from design altitude + Mach.
    # This is the same logic as set_ambient_conditions! in size_aircraft.jl,
    # reproduced here so the harness does not depend on that unexported helper.
    # -----------------------------------------------------------------------
    ΔTatmos = ac.parm[imDeltaTatm]
    alt_m   = ac.para[iaalt, ip, imission]
    as      = atmos(alt_m, ΔTatmos)
    Mach    = ac.para[iaMach, ip, imission]

    ac.pare[iep0,   ip, imission] = as.p
    ac.pare[ieT0,   ip, imission] = as.T
    ac.pare[iea0,   ip, imission] = as.a
    ac.pare[ierho0, ip, imission] = as.ρ
    ac.pare[iemu0,  ip, imission] = as.μ
    ac.pare[ieM0,   ip, imission] = Mach
    ac.pare[ieu0,   ip, imission] = Mach * as.a
    ac.para[iaReunit, ip, imission] = Mach * as.a * as.ρ / as.μ

    # -----------------------------------------------------------------------
    # Run engine design-point sizing (no aircraft weight / aero loop)
    # -----------------------------------------------------------------------
    tfwrap!(ac, "design", imission, ip, true)

    # -----------------------------------------------------------------------
    # Read converged pare column → typed EngineState
    # -----------------------------------------------------------------------
    eng = EngineState{Float64}()
    pare_to_engine_state!(eng, view(ac.pare, :, ip, imission))
    return eng
end
