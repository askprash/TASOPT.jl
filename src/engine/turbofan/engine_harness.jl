"""
engine_harness.jl — Engine-standalone runners (single-point and multi-point sweep).

Provides:
- `pare_to_engine_state!`: public mapping from a `pare` column to a typed
  `EngineState`.
- `run_engine_design_point`: single-point design sizing runner — accepts an
  aircraft object, initialises ambient conditions, runs `tfsize!` via
  `tfwrap!`, and returns an `EngineState`.
- `run_engine_sweep`: multi-point off-design runner — iterates the engine
  across a range of mission points, writes results into `ac.missions[imission]`,
  and returns that `Mission{Float64}`.
- `write_sweep_csv`: serialise a `Mission` sweep to CSV.
- `write_sweep_toml`: serialise a `Mission` sweep to TOML regression baseline.
- `SweepResult`: **deprecated** tabular container (use `Mission{T}` instead).
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
- **Cooling design state**: `design.epsrow` (per-row cooling bypass ratios),
  `design.Tmrow` (blade metal temperatures), `design.fc` (total cooling fraction).

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

    # -----------------------------------------------------------------------
    # Cooling design state — always populated so tfcalc! can read from typed
    # state rather than indexing pare[ieepsc1..4] directly.
    # epsrow[i] = ṁ_cool,i / ṁ_core (cooling bypass ratio for blade row i)
    # Tmrow[i]  = blade metal temperature [K] for blade row i
    # fc        = total cooling mass-flow fraction (normalised by design core)
    # -----------------------------------------------------------------------
    eng.design.epsrow = SVector{4,Float64}(
        pare[ieepsc1], pare[ieepsc2], pare[ieepsc3], pare[ieepsc4])
    eng.design.Tmrow  = SVector{4,Float64}(
        pare[ieTmet1], pare[ieTmet2], pare[ieTmet3], pare[ieTmet4])
    eng.design.fc     = pare[iefc]

    # -----------------------------------------------------------------------
    # Design map anchors — design-point PRs, corrected flows, corrected speeds
    # (tasopt-j9l.23): populated so off-design tfcalc! reads from typed state
    # rather than indexing pare[iepifD..iepiltD] / pare[iembfD..iembltD] /
    # pare[ieNbfD..ieNbltD] directly.
    # -----------------------------------------------------------------------
    eng.design.pifD  = pare[iepifD]
    eng.design.pilcD = pare[iepilcD]
    eng.design.pihcD = pare[iepihcD]
    eng.design.pihtD = pare[iepihtD]
    eng.design.piltD = pare[iepiltD]

    eng.design.mbfD  = pare[iembfD]
    eng.design.mblcD = pare[iemblcD]
    eng.design.mbhcD = pare[iembhcD]
    eng.design.mbhtD = pare[iembhtD]
    eng.design.mbltD = pare[iembltD]

    eng.design.NbfD  = pare[ieNbfD]
    eng.design.NblcD = pare[ieNblcD]
    eng.design.NbhcD = pare[ieNbhcD]
    eng.design.NbhtD = pare[ieNbhtD]
    eng.design.NbltD = pare[ieNbltD]

    # -----------------------------------------------------------------------
    # Design component flow areas (tasopt-j9l.53)
    # These frozen sizing-point areas are read so off-design tfcalc! can
    # consume them from typed DesignState rather than indexing pare directly.
    # -----------------------------------------------------------------------
    eng.design.A2  = pare[ieA2]
    eng.design.A25 = pare[ieA25]
    eng.design.A5  = pare[ieA5]
    eng.design.A7  = pare[ieA7]

    # -----------------------------------------------------------------------
    # Performance rollup scalars (tasopt-j9l.52)
    # -----------------------------------------------------------------------
    eng.TSFC  = pare[ieTSFC]
    eng.Fe    = pare[ieFe]
    eng.Fsp   = pare[ieFsp]
    eng.BPR   = pare[ieBPR]
    eng.mfuel = pare[iemfuel]

    # -----------------------------------------------------------------------
    # Component adiabatic efficiencies (tasopt-j9l.63.1)
    # -----------------------------------------------------------------------
    eng.etaf  = pare[ieetaf]
    eng.etalc = pare[ieetalc]
    eng.etahc = pare[ieetahc]
    eng.etaht = pare[ieetaht]
    eng.etalt = pare[ieetalt]

    # -----------------------------------------------------------------------
    # Overall propulsion efficiencies (tasopt-j9l.63.2)
    # No pare backing — derived from scalar pare fields that ARE stored.
    # Formula: eta_overall = gee * u0 / (TSFC * hfuel)   [≡ Fe*u0/(mdotf*hfuel)]
    #          eta_prop    = 2 / (2 + Fsp)               [mixed-jet momentum]
    #          eta_thermal = eta_overall / eta_prop
    # Guard: zero at ground-idle where TSFC ≤ 0 or Fsp ≤ 0.
    # -----------------------------------------------------------------------
    let _TSFC  = pare[ieTSFC],
        _Fsp   = pare[ieFsp],
        _u0    = pare[ieu0],
        _hfuel = pare[iehfuel]
        if _TSFC > 0 && _Fsp > 0
            _eta_ov      = gee * _u0 / (_TSFC * _hfuel)
            _eta_pr      = 2 / (2 + _Fsp)
            eng.eta_overall = _eta_ov
            eng.eta_prop    = _eta_pr
            eng.eta_thermal = _eta_ov / _eta_pr
        else
            eng.eta_overall = zero(eltype(pare))
            eng.eta_prop    = zero(eltype(pare))
            eng.eta_thermal = zero(eltype(pare))
        end
    end

    return eng
end

# ---------------------------------------------------------------------------
# engine_state_to_pare! — inverse of pare_to_engine_state!
# ---------------------------------------------------------------------------

"""
    engine_state_to_pare!(eng, pare) -> eng

Write the thermodynamic station state of an `EngineState` back into a single
column of the `pare` array.  This is the inverse of [`pare_to_engine_state!`](@ref).

Only fields that the legacy `pare` array actually stores are written (same
coverage as `pare_to_engine_state!`):

- **Ambient scalars**: `M0`, `T0`, `p0`, `a0`.
- **Total state** (Tt, ht, pt, cpt, Rt): stations 0, 18, 19, 2, 21, 25, 3,
  4 (ht/pt/cpt/Rt only — Tt4 is an input, not written), 41, 45, 49, 5, 7.
- **Static state** (ps, Ts, Rs, cps, u): stations 2, 25, 5, 6, 7, 8.
- **Area** (A): stations 2, 25, 5, 6, 7, 8, 9.
- **Mass flow** (mdot): station 2 (`mcore`).
- **Station 9**: Tt, pt, u, A.
- **Performance rollup**: `TSFC`, `Fsp`, `BPR`, `mfuel`.  `Fe` is NOT written here
  (it is mode-dependent: a conditional bare write handles FixedTt4OffDes output).

Fields not in `pare` (stations 19c, 25c, 4a, 49c; `hs`, `ss`) are not written.
"""
function engine_state_to_pare!(eng::EngineState, pare)
    # Ambient scalars
    pare[ieM0] = eng.M0
    pare[ieT0] = eng.T0
    pare[iep0] = eng.p0
    pare[iea0] = eng.a0

    # Station 0 — freestream
    _write_total!(pare, ieTt0, ieht0, iept0, iecpt0, ieRt0, eng.st0)
    pare[ieu0] = eng.st0.u

    # Station 18 — FanFaceOuter
    _write_total!(pare, ieTt18, ieht18, iept18, iecpt18, ieRt18, eng.st18)

    # Station 19 — FanFaceLPC
    _write_total!(pare, ieTt19, ieht19, iept19, iecpt19, ieRt19, eng.st19)

    # Station 2 — FanFaceFan  (total + static + area + mdot)
    _write_total!(pare, ieTt2, ieht2, iept2, iecpt2, ieRt2, eng.st2)
    _write_static!(pare, iep2, ieT2, ieR2, iecp2, ieu2, eng.st2)
    pare[ieA2]     = eng.st2.A
    pare[iemcore]  = eng.st2.mdot

    # Station 21 — FanExit
    _write_total!(pare, ieTt21, ieht21, iept21, iecpt21, ieRt21, eng.st21)

    # Station 25 — LPCExit  (total + static + area)
    _write_total!(pare, ieTt25, ieht25, iept25, iecpt25, ieRt25, eng.st25)
    _write_static!(pare, iep25, ieT25, ieR25, iecp25, ieu25, eng.st25)
    pare[ieA25] = eng.st25.A

    # Station 3 — HPCExit
    _write_total!(pare, ieTt3, ieht3, iept3, iecpt3, ieRt3, eng.st3)

    # Station 4 — CombustorExit
    # NOTE: Tt4 is an INPUT to tfsize! (not written back); only derived
    # quantities ht4/pt4/cpt4/Rt4 are outputs.
    pare[ieht4]  = eng.st4.ht
    pare[iept4]  = eng.st4.pt
    pare[iecpt4] = eng.st4.cpt
    pare[ieRt4]  = eng.st4.Rt

    # Station 41 — TurbineInlet
    _write_total!(pare, ieTt41, ieht41, iept41, iecpt41, ieRt41, eng.st41)

    # Station 45 — HPTExit
    _write_total!(pare, ieTt45, ieht45, iept45, iecpt45, ieRt45, eng.st45)

    # Station 49 — LPTExit
    _write_total!(pare, ieTt49, ieht49, iept49, iecpt49, ieRt49, eng.st49)

    # Station 5 — CoreNozzle  (total + static + area)
    _write_total!(pare, ieTt5, ieht5, iept5, iecpt5, ieRt5, eng.st5)
    _write_static!(pare, iep5, ieT5, ieR5, iecp5, ieu5, eng.st5)
    pare[ieA5] = eng.st5.A

    # Station 6 — CoreNozzleExit  (static + area; no total in pare)
    _write_static!(pare, iep6, ieT6, ieR6, iecp6, ieu6, eng.st6)
    pare[ieA6] = eng.st6.A

    # Station 7 — FanNozzle  (total + static + area)
    _write_total!(pare, ieTt7, ieht7, iept7, iecpt7, ieRt7, eng.st7)
    _write_static!(pare, iep7, ieT7, ieR7, iecp7, ieu7, eng.st7)
    pare[ieA7] = eng.st7.A

    # Station 8 — FanNozzleExit  (static + area; no total in pare)
    _write_static!(pare, iep8, ieT8, ieR8, iecp8, ieu8, eng.st8)
    pare[ieA8] = eng.st8.A

    # Station 9 — OfftakeDisch
    pare[ieTt9] = eng.st9.Tt
    pare[iept9] = eng.st9.pt
    pare[ieu9]  = eng.st9.u
    pare[ieA9]  = eng.st9.A

    # Cooling — per-row cooling fractions (ieepsc1..4) and blade metal
    # temperatures (ieTmet1..4).  Both arrays are always populated in the EXIT
    # blocks (sizing lines 376-378; off-design just before this call), so we
    # can write them unconditionally here.
    # NOTE: pare[iefc] is intentionally NOT written here.  In FixedTmetal mode
    # it is written explicitly in the off-design EXIT block (and sizing EXIT
    # block) because fc is a computed output only in that mode.  In
    # FixedCoolingFlowRatio mode fc is not recomputed and pare[iefc] retains
    # the initialization value — matching historical behaviour.
    for icrow in 1:length(eng.design.epsrow)
        pare[ieepsc1+icrow-1] = eng.design.epsrow[icrow]
        pare[ieTmet1+icrow-1] = eng.design.Tmrow[icrow]
    end

    # Performance rollup scalars (tasopt-j9l.52)
    # NOTE: pare[ieFe] is intentionally NOT written here.
    # Fe is a MODE-DEPENDENT input/output:
    #   - CalcMode.Sizing:          Fe = INPUT (design thrust target)
    #   - CalcMode.FixedTt4OffDes:  Fe = OUTPUT (computed thrust)
    #   - CalcMode.FixedFeOffDes:   Fe = INPUT (required thrust from mission optimization)
    # Writing Fe unconditionally would corrupt the required thrust in FixedFeOffDes
    # mode (Newton convergence gives Fe_returned ≠ Fe_required at floating-point level).
    # The conditional bare write `if FixedTt4OffDes: pare[ieFe] = Fe` handles the
    # off-design output case, exactly mirroring the ieTt4 pattern for FixedFeOffDes.
    pare[ieTSFC]  = eng.TSFC
    pare[ieFsp]   = eng.Fsp
    pare[ieBPR]   = eng.BPR
    pare[iemfuel] = eng.mfuel

    # Component adiabatic efficiencies (tasopt-j9l.63.1)
    pare[ieetaf]  = eng.etaf
    pare[ieetalc] = eng.etalc
    pare[ieetahc] = eng.etahc
    pare[ieetaht] = eng.etaht
    pare[ieetalt] = eng.etalt

    # Overall propulsion efficiencies (tasopt-j9l.63.2)
    # No pare backing — intentionally not written to pare.

    return eng
end

"""
    _write_total!(pare, iTt, iht, ipt, icpt, iRt, fs)

Write the five total-state scalars from `fs` into a `pare` slice.
Inverse of `_fill_total!`.
"""
@inline function _write_total!(pare, iTt::Int, iht::Int, ipt::Int, icpt::Int, iRt::Int,
                                fs::FlowStation)
    pare[iTt]  = fs.Tt
    pare[iht]  = fs.ht
    pare[ipt]  = fs.pt
    pare[icpt] = fs.cpt
    pare[iRt]  = fs.Rt
end

"""
    _write_static!(pare, ip, iT, iR, icp, iu, fs)

Write the five static-state scalars (p, T, R, cp, u) from `fs` into a `pare`
slice.  Inverse of `_fill_static!`.
"""
@inline function _write_static!(pare, ip::Int, iT::Int, iR::Int, icp::Int, iu::Int,
                                 fs::FlowStation)
    pare[ip]  = fs.ps
    pare[iT]  = fs.Ts
    pare[iR]  = fs.Rs
    pare[icp] = fs.cps
    pare[iu]  = fs.u
end

# ---------------------------------------------------------------------------
# design_state_to_pare! — write DesignState scalars back to pare
# ---------------------------------------------------------------------------

"""
    design_state_to_pare!(ds, pare) -> ds

Write the frozen design-point scalars from a `DesignState` back into the
`pare` slice.  These are the quantities set once during `tfsize!` (design
pressure ratios, corrected flows, corrected spool speeds, flow areas) and
subsequently read by `tfoper!` for every off-design point.

The cooling arrays (`epsrow`, `Tmrow`, `fc`) are written by
`engine_state_to_pare!` via the `EngineState.design` fields, not here.
"""
function design_state_to_pare!(ds::DesignState, pare)
    # Flow areas
    pare[ieA2]  = ds.A2
    pare[ieA25] = ds.A25
    pare[ieA5]  = ds.A5
    pare[ieA7]  = ds.A7

    # Design-point corrected spool speeds
    pare[ieNbfD]  = ds.NbfD
    pare[ieNblcD] = ds.NblcD
    pare[ieNbhcD] = ds.NbhcD
    pare[ieNbhtD] = ds.NbhtD
    pare[ieNbltD] = ds.NbltD

    # Design-point corrected mass flows
    pare[iembfD]  = ds.mbfD
    pare[iemblcD] = ds.mblcD
    pare[iembhcD] = ds.mbhcD
    pare[iembhtD] = ds.mbhtD
    pare[iembltD] = ds.mbltD

    # Design-point pressure ratios
    pare[iepifD]  = ds.pifD
    pare[iepilcD] = ds.pilcD
    pare[iepihcD] = ds.pihcD
    pare[iepihtD] = ds.pihtD
    pare[iepiltD] = ds.piltD

    return ds
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

# ---------------------------------------------------------------------------
# SweepResult
# ---------------------------------------------------------------------------

"""
    SweepResult{T<:AbstractFloat}

!!! warning "Deprecated"
    `SweepResult` is deprecated.  [`run_engine_sweep`](@ref) now returns a
    `Mission{T}` directly and writes results into `ac.missions[imission]`.
    Use `mission.points[ip].engine` to access per-point engine state.
    `SweepResult` will be removed in a future release.

Tabular container for the output of [`run_engine_sweep`](@ref).  Holds one
[`EngineState`](@ref) per mission point together with the key engine
performance scalars extracted from `pare`.

## Fields

| Field       | Size  | Unit    | Description                              |
|:------------|:------|:--------|:-----------------------------------------|
| `ip_indices`| `n`   | —       | Mission-point indices (`ip` values)       |
| `ip_labels` | `n`   | —       | Human-readable labels (e.g. `"cruise1"`) |
| `alt`       | `n`   | m       | Altitude at each point                   |
| `Mach`      | `n`   | —       | Mach number at each point                |
| `engines`   | `n`   | —       | `EngineState` per point                  |
| `Fe`        | `n`   | N       | Net thrust per engine                    |
| `TSFC`      | `n`   | kg/N/s  | Thrust-specific fuel consumption         |
| `BPR`       | `n`   | —       | Bypass ratio                             |
| `mcore`     | `n`   | kg/s    | Core mass flow rate (per engine)         |
| `mdotf`     | `n`   | kg/s    | Fuel mass flow rate (total, all engines; from `pare[iemfuel]`) |
"""
struct SweepResult{T<:AbstractFloat}
    ip_indices ::Vector{Int}
    ip_labels  ::Vector{String}
    alt        ::Vector{T}
    Mach       ::Vector{T}
    engines    ::Vector{EngineState{T}}
    Fe         ::Vector{T}
    TSFC       ::Vector{T}
    BPR        ::Vector{T}
    mcore      ::Vector{T}
    mdotf      ::Vector{T}
end

# ---------------------------------------------------------------------------
# run_engine_sweep
# ---------------------------------------------------------------------------

"""
    run_engine_sweep(ac; imission=1, ip_range=ipstatic:ipdescentn,
                    initializes_engine=false) -> Mission{Float64}

Run the turbofan off-design performance routine across every mission point in
`ip_range`, write the converged engine state into `ac.missions[imission]`, and
return that `Mission{Float64}` directly.

The aircraft must already be sized with `size_aircraft!` (or equivalent)
before calling this function.  In particular:
- `pare[ieFe, ip, imission]` must hold the required thrust at each mission
  point (set by `size_aircraft!` ↔ `fly_mission!`).
- Ambient conditions in `pare` (T0, p0, a0, M0, …) must be populated.
- `ac.missions[imission]` must be pre-allocated (done by `read_input.jl`).

For climb and ground-roll points (`ipstatic:ipclimbn`) the engine is operated
at maximum turbine-inlet temperature (`CalcMode.FixedTt4OffDes`); for all
other points it operates at the thrust target (`CalcMode.FixedFeOffDes`).

## Arguments
- `ac`: a sized `aircraft` object.
- `imission::Int`: mission index (default `1`).
- `ip_range`: iterable of mission-point indices (default `ipstatic:ipdescentn`,
  i.e., all 16 regular mission points).
- `initializes_engine::Bool`: passed to `tfwrap!` for each point.
  Default `false` — use the converged `pare` state as the initial guess.
  Pass `true` to reinitialise the Newton iteration from scratch at every point
  (slower, but independent of the current `pare` state).

## Returns
`ac.missions[imission]::Mission{Float64}` — the mission whose `points[ip].engine`
fields have been populated for every `ip` in `ip_range`.

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
size_aircraft!(ac; printiter=false)
mission = TASOPT.engine.run_engine_sweep(ac)
TASOPT.engine.write_sweep_csv("engine_sweep.csv", mission, ipstatic:ipdescentn, ac)
```
"""
function run_engine_sweep(ac;
                          imission::Int        = 1,
                          ip_range             = ipstatic:ipdescentn,
                          initializes_engine::Bool = false)

    mission = ac.missions[imission]

    for ip in ip_range
        # Run engine off-design at this mission point
        tfwrap!(ac, "off_design", imission, ip, initializes_engine)

        # Map converged pare column → typed EngineState in mission
        pare_to_engine_state!(mission.points[ip].engine,
                               view(ac.pare, :, ip, imission))
    end

    return mission
end

# ---------------------------------------------------------------------------
# write_sweep_csv
# ---------------------------------------------------------------------------

const _SWEEP_CSV_STATIONS = (
    ("0",  :st0),  ("2",  :st2),  ("3",  :st3),
    ("4",  :st4),  ("41", :st41), ("45", :st45),
    ("49", :st49), ("5",  :st5),  ("7",  :st7))

"""
    write_sweep_csv(io::IO, mission::Mission, ip_range, ac; imission=1)
    write_sweep_csv(path::AbstractString, mission::Mission, ip_range, ac; imission=1)

Write a [`Mission`](@ref) sweep to CSV format for the mission points in `ip_range`.

Each row corresponds to one mission point.  Columns:
- Metadata: `ip`, `label`, `alt_m`, `Mach`
- Engine performance: `Fe_N`, `TSFC_kg_Ns`, `BPR`, `mcore_kg_s`, `mfuel_kg_s`
- Station totals at key stations (0, 2, 3, 4, 41, 45, 49, 5, 7):
  `Tt<N>_K`, `pt<N>_Pa`, `ht<N>_J_kg`
- Station exit velocities at nozzle stations 5 and 7:
  `u5_m_s`, `u7_m_s`

The `path` overload opens the file, writes, and closes it.
"""
function write_sweep_csv(io::IO, mission, ip_range, ac; imission::Int=1)
    # ----- Header -----
    header_parts = [
        "ip", "label", "alt_m", "Mach",
        "Fe_N", "TSFC_kg_Ns", "BPR", "mcore_kg_s", "mfuel_kg_s",
    ]
    for (sname, _) in _SWEEP_CSV_STATIONS
        push!(header_parts, "Tt$(sname)_K", "pt$(sname)_Pa", "ht$(sname)_J_kg")
    end
    push!(header_parts, "u5_m_s", "u7_m_s")
    println(io, join(header_parts, ","))

    # ----- Rows -----
    for ip in ip_range
        eng = mission.points[ip].engine
        lbl = (ip <= length(ip_labels)) ? ip_labels[ip] : string(ip)
        row = [
            string(ip),
            lbl,
            @sprintf("%.2f", ac.para[iaalt,  ip, imission]),
            @sprintf("%.6g", ac.para[iaMach, ip, imission]),
            @sprintf("%.6g", eng.Fe),
            @sprintf("%.6g", eng.TSFC),
            @sprintf("%.6g", eng.BPR),
            @sprintf("%.6g", eng.st2.mdot),
            @sprintf("%.6g", eng.mfuel),
        ]
        for (_, stfield) in _SWEEP_CSV_STATIONS
            fs = getfield(eng, stfield)
            push!(row,
                  @sprintf("%.6g", fs.Tt),
                  @sprintf("%.6g", fs.pt),
                  @sprintf("%.6g", fs.ht))
        end
        push!(row,
              @sprintf("%.6g", eng.st5.u),
              @sprintf("%.6g", eng.st7.u))
        println(io, join(row, ","))
    end
    return nothing
end

function write_sweep_csv(path::AbstractString, mission, ip_range, ac; imission::Int=1)
    open(path, "w") do io
        write_sweep_csv(io, mission, ip_range, ac; imission=imission)
    end
    return nothing
end

"""
    write_sweep_csv(io::IO, result::SweepResult)
    write_sweep_csv(path::AbstractString, result::SweepResult)

!!! warning "Deprecated"
    This overload is deprecated.  Use `write_sweep_csv(io, mission, ip_range, ac)`
    with the `Mission{T}` returned by [`run_engine_sweep`](@ref).

Write a [`SweepResult`](@ref) to CSV format.
"""
function write_sweep_csv(io::IO, result::SweepResult)
    # ----- Header -----
    header_parts = [
        "ip", "label", "alt_m", "Mach",
        "Fe_N", "TSFC_kg_Ns", "BPR", "mcore_kg_s", "mfuel_kg_s",
    ]
    for (sname, _) in _SWEEP_CSV_STATIONS
        push!(header_parts, "Tt$(sname)_K", "pt$(sname)_Pa", "ht$(sname)_J_kg")
    end
    push!(header_parts, "u5_m_s", "u7_m_s")
    println(io, join(header_parts, ","))

    # ----- Rows -----
    n = length(result.ip_indices)
    for k in 1:n
        eng = result.engines[k]
        row = [
            string(result.ip_indices[k]),
            result.ip_labels[k],
            @sprintf("%.2f", result.alt[k]),
            @sprintf("%.6g", result.Mach[k]),
            @sprintf("%.6g", result.Fe[k]),
            @sprintf("%.6g", result.TSFC[k]),
            @sprintf("%.6g", result.BPR[k]),
            @sprintf("%.6g", result.mcore[k]),
            @sprintf("%.6g", result.mdotf[k]),
        ]
        for (_, stfield) in _SWEEP_CSV_STATIONS
            fs = getfield(eng, stfield)
            push!(row,
                  @sprintf("%.6g", fs.Tt),
                  @sprintf("%.6g", fs.pt),
                  @sprintf("%.6g", fs.ht))
        end
        push!(row,
              @sprintf("%.6g", eng.st5.u),
              @sprintf("%.6g", eng.st7.u))
        println(io, join(row, ","))
    end
    return nothing
end

function write_sweep_csv(path::AbstractString, result::SweepResult)
    open(path, "w") do io
        write_sweep_csv(io, result)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# write_sweep_toml — diffable full-fidelity station baseline
# ---------------------------------------------------------------------------

import TOML

# Stations serialised in physical flow-path order (mirrors _STATION_DUMP_ORDER).
const _TOML_STATION_ORDER = (
    ("0",   "Freestream",     :st0),
    ("2",   "FanFaceFan",     :st2),
    ("18",  "FanFaceOuter",   :st18),
    ("19",  "FanFaceLPC",     :st19),
    ("19c", "PreCoolerOut",   :st19c),
    ("21",  "FanExit",        :st21),
    ("25",  "LPCExit",        :st25),
    ("25c", "InterCoolerOut", :st25c),
    ("3",   "HPCExit",        :st3),
    ("4",   "CombustorExit",  :st4),
    ("4a",  "CoolMixInlet",   :st4a),
    ("41",  "TurbineInlet",   :st41),
    ("45",  "HPTExit",        :st45),
    ("49",  "LPTExit",        :st49),
    ("49c", "RegenCoolerOut", :st49c),
    ("5",   "CoreNozzle",     :st5),
    ("6",   "CoreNozzleExit", :st6),
    ("7",   "FanNozzle",      :st7),
    ("8",   "FanNozzleExit",  :st8),
    ("9",   "OfftakeDisch",   :st9),
)

"""
    _station_to_dict(fs::FlowStation) -> Dict{String,Float64}

Serialise the populated scalar fields of a `FlowStation` to a plain
`Dict{String,Float64}`.  Fields not written by `pare_to_engine_state!`
(hs, ss, st, alpha) are excluded to keep the baseline compact and to
avoid noisy diffs on zero-valued entries.
"""
function _station_to_dict(fs::FlowStation)
    Dict{String,Float64}(
        "Tt"   => Float64(fs.Tt),
        "ht"   => Float64(fs.ht),
        "pt"   => Float64(fs.pt),
        "cpt"  => Float64(fs.cpt),
        "Rt"   => Float64(fs.Rt),
        "Ts"   => Float64(fs.Ts),
        "ps"   => Float64(fs.ps),
        "cps"  => Float64(fs.cps),
        "Rs"   => Float64(fs.Rs),
        "u"    => Float64(fs.u),
        "A"    => Float64(fs.A),
        "mdot" => Float64(fs.mdot),
    )
end

"""
    write_sweep_toml(io::IO, mission::Mission, ip_range, ac; imission=1)
    write_sweep_toml(path::AbstractString, mission::Mission, ip_range, ac; imission=1)

Write a [`Mission`](@ref) sweep to TOML format suitable for use as a
regression baseline fixture.

Each `[[points]]` table contains:
- Metadata: `ip`, `label`, `alt_m`, `Mach`.
- Engine performance: `Fe_N`, `TSFC_kg_Ns`, `BPR`, `mcore_kg_s`,
  `mfuel_kg_s`, `M0`, `T0_K`, `p0_Pa`, `a0_m_s`.
- Station subtables (e.g. `[points.stations.st0]`) with twelve scalar
  fields: `Tt`, `ht`, `pt`, `cpt`, `Rt`, `Ts`, `ps`, `cps`, `Rs`,
  `u`, `A`, `mdot`.

The TOML format is line-oriented and deterministic, making it suitable
for `git diff` review when baseline values change.
"""
function write_sweep_toml(io::IO, mission, ip_range, ac; imission::Int=1)
    ips = collect(Int, ip_range)
    n = length(ips)
    points = Vector{Dict{String,Any}}(undef, n)
    for (k, ip) in enumerate(ips)
        eng = mission.points[ip].engine
        lbl = (ip <= length(ip_labels)) ? ip_labels[ip] : string(ip)
        stations = Dict{String,Any}()
        for (_, _, stfld) in _TOML_STATION_ORDER
            stations[String(stfld)] = _station_to_dict(getfield(eng, stfld))
        end
        points[k] = Dict{String,Any}(
            "ip"          => ip,
            "label"       => lbl,
            "alt_m"       => Float64(ac.para[iaalt,  ip, imission]),
            "Mach"        => Float64(ac.para[iaMach, ip, imission]),
            "Fe_N"        => Float64(eng.Fe),
            "TSFC_kg_Ns"  => Float64(eng.TSFC),
            "BPR"         => Float64(eng.BPR),
            "mcore_kg_s"  => Float64(eng.st2.mdot),
            "mfuel_kg_s"  => Float64(eng.mfuel),
            "M0"          => Float64(eng.M0),
            "T0_K"        => Float64(eng.T0),
            "p0_Pa"       => Float64(eng.p0),
            "a0_m_s"      => Float64(eng.a0),
            "stations"    => stations,
        )
    end
    data = Dict{String,Any}(
        "metadata" => Dict{String,Any}(
            "n_points" => n,
            "aircraft" => "default",
            "description" => "Engine sweep baseline — default TASOPT aircraft.",
        ),
        "points" => points,
    )
    TOML.print(io, data)
    return nothing
end

function write_sweep_toml(path::AbstractString, mission, ip_range, ac; imission::Int=1)
    open(path, "w") do io
        write_sweep_toml(io, mission, ip_range, ac; imission=imission)
    end
    return nothing
end

"""
    write_sweep_toml(io::IO, result::SweepResult)
    write_sweep_toml(path::AbstractString, result::SweepResult)

!!! warning "Deprecated"
    This overload is deprecated.  Use `write_sweep_toml(io, mission, ip_range, ac)`
    with the `Mission{T}` returned by [`run_engine_sweep`](@ref).

Write a [`SweepResult`](@ref) to TOML format.
"""
function write_sweep_toml(io::IO, result::SweepResult)
    n = length(result.ip_indices)
    points = Vector{Dict{String,Any}}(undef, n)
    for k in 1:n
        eng = result.engines[k]
        stations = Dict{String,Any}()
        for (_, _, stfld) in _TOML_STATION_ORDER
            stations[String(stfld)] = _station_to_dict(getfield(eng, stfld))
        end
        points[k] = Dict{String,Any}(
            "ip"          => result.ip_indices[k],
            "label"       => result.ip_labels[k],
            "alt_m"       => Float64(result.alt[k]),
            "Mach"        => Float64(result.Mach[k]),
            "Fe_N"        => Float64(result.Fe[k]),
            "TSFC_kg_Ns"  => Float64(result.TSFC[k]),
            "BPR"         => Float64(result.BPR[k]),
            "mcore_kg_s"  => Float64(result.mcore[k]),
            "mfuel_kg_s"  => Float64(result.mdotf[k]),
            "M0"          => Float64(eng.M0),
            "T0_K"        => Float64(eng.T0),
            "p0_Pa"       => Float64(eng.p0),
            "a0_m_s"      => Float64(eng.a0),
            "stations"    => stations,
        )
    end
    data = Dict{String,Any}(
        "metadata" => Dict{String,Any}(
            "n_points" => n,
            "aircraft" => "default",
            "description" => "Engine sweep baseline — default TASOPT aircraft.",
        ),
        "points" => points,
    )
    TOML.print(io, data)
    return nothing
end

function write_sweep_toml(path::AbstractString, result::SweepResult)
    open(path, "w") do io
        write_sweep_toml(io, result)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# regenerate_engine_baseline
# ---------------------------------------------------------------------------

"""
    const ENGINE_BASELINE_PATH

Default path for the engine-sweep regression baseline fixture.
Points to `test/fixtures/engine_sweep_baseline.toml` relative to the
repository root.
"""
const ENGINE_BASELINE_PATH = joinpath(
    @__DIR__, "..", "..", "..", "test", "fixtures", "engine_sweep_baseline.toml")

"""
    regenerate_engine_baseline(ac; path=ENGINE_BASELINE_PATH)

**Auditable action** — regenerates the engine-sweep regression baseline
from a sized aircraft object and writes it to `path`.

Run the full off-design sweep on `ac` (which must already be sized with
`size_aircraft!`) and serialise the result to TOML via `write_sweep_toml`.

## ⚠  Warning

Regenerating the baseline changes the numerical reference used by the
engine regression test (`tasopt-j9l.5`).  Only call this function when a
physics or model change is intentional and the new values have been
verified.  Commit the updated baseline with an explanatory git message —
reviewers will see the numerical diff in the pull request.

## Usage
```julia
using TASOPT
ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)

# Regenerate to the default test-fixture path
TASOPT.engine.regenerate_engine_baseline(ac)

# Or write to a custom path for inspection first
TASOPT.engine.regenerate_engine_baseline(ac; path="/tmp/my_baseline.toml")
```
"""
function regenerate_engine_baseline(ac; path::AbstractString=ENGINE_BASELINE_PATH)
    @warn """
    *** REGENERATING ENGINE SWEEP BASELINE ***
    Output path: $(abspath(path))

    This rewrites the numerical reference for engine regression tests.
    Only proceed when the change is intentional and has been verified.
    Commit the updated baseline file with an explanatory git message so
    reviewers can audit the numerical diff.
    """
    mkpath(dirname(abspath(path)))
    mission = run_engine_sweep(ac)
    write_sweep_toml(path, mission, ipstatic:ipdescentn, ac)
    n = ipdescentn - ipstatic + 1
    @info "Baseline written: $(abspath(path))  ($n mission points)"
    return abspath(path)
end
