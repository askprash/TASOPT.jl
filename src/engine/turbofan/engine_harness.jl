"""
engine_harness.jl — Engine-standalone runners (single-point and multi-point sweep).

Provides:
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
# (still used by ducted_fan_harness.jl)
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
# sync_freestream_to_pare!
# ---------------------------------------------------------------------------

"""
    sync_freestream_to_pare!(eng, pare)

Write the seven ambient freestream scalars (`p0`, `T0`, `a0`, `rho0`, `mu0`,
`M0`, `u0`) from an `EngineState` into the corresponding `pare` column entries.

Used by mission_iteration.jl for mission points that do not go through a full
engine calculation (iptakeoff, ipcutback, intermediate cruise), so that
`pare[iep0..ieu0, ip]` stays consistent with typed state between iterations.
"""
function sync_freestream_to_pare!(eng::EngineState, pare)
    pare[iep0]   = eng.p0
    pare[ieT0]   = eng.T0
    pare[iea0]   = eng.a0
    pare[ierho0] = eng.rho0
    pare[iemu0]  = eng.mu0
    pare[ieM0]   = eng.M0
    pare[ieu0]   = eng.st0.u
end

# ---------------------------------------------------------------------------
# sync_design_scalars_to_pare! — write compressor-map anchor scalars to pare
# ---------------------------------------------------------------------------

"""
    sync_design_scalars_to_pare!(ds, pare) -> ds

Write the compressor-map anchor scalars from a `DesignState` to a `pare`
slice: flow areas (A2, A25, A5, A7), corrected spool speeds (NbfD..NbltD),
corrected mass flows (mbfD..mbltD), and pressure ratios (pifD..piltD).

Unlike `design_state_to_pare!` (deleted), this does **not** write `ruc` or `M4a`.
Those cooling mixing constants are initialised from the design-point bare-pare
write (`sync_design_scalars_to_pare!` via `tfcalc!`) and must not be overwritten
for non-design mission points during the per-point broadcast in `tfwrap!`
(tasopt-j9l.45.14.2).
"""
function sync_design_scalars_to_pare!(ds::DesignState, pare)
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
# sync_cooling_scalars_to_pare! — write cooling design scalars to pare
# ---------------------------------------------------------------------------

"""
    sync_cooling_scalars_to_pare!(ds, pare) -> ds

Write the blade-row cooling bypass fractions (`epsrow`) and total cooling
mass-flow fraction (`fc`) from a `DesignState` to a `pare` slice.

Called after `FixedTmetal` cooling sizing to keep bare pare in sync with the
typed EngineState (tasopt-j9l.45.14.2).
"""
function sync_cooling_scalars_to_pare!(ds::DesignState, pare)
    for icrow in 1:length(ds.epsrow)
        pare[ieepsc1+icrow-1] = ds.epsrow[icrow]
    end
    pare[iefc] = ds.fc
    return ds
end

# ---------------------------------------------------------------------------
# engine_state_to_pare_vec — read typed EngineState → bare-pare-layout vector
# ---------------------------------------------------------------------------

"""
    engine_state_to_pare_vec(eng::EngineState) -> Vector{Float64}

Return a `ietotal`-length `Float64` vector whose layout mirrors the bare-pare
`ie*` index table (index.inc).  Each slot that has a typed-state equivalent is
filled from `eng` or `eng.design`; slots with no typed-state mapping (NOx
indices, radiator ε, ieDi, etc.) are set to `NaN`.

Used by `output_csv` to eliminate bare-pare reads from IO without requiring
a parallel `ac.pare` access.
"""
function engine_state_to_pare_vec(eng::EngineState)::Vector{Float64}
    ds = eng.design
    v  = fill(NaN, ietotal)

    # ---- scalars: per-point operating outputs ----
    v[iehfuel]  = eng.hfuel
    v[ieTfuel]  = eng.Tfuel
    v[ieff]     = eng.ff
    # ---- frozen design inputs (iepid..iepitn) ----
    v[iepid]    = ds.pid
    v[iepib]    = ds.pib
    v[iepifn]   = ds.pifn
    v[iepitn]   = ds.pitn
    # ---- performance output ----
    v[ieBPR]    = eng.BPR
    # ---- frozen polytropic efficiencies ----
    v[ieepolf]  = ds.epolf
    v[ieepollc] = ds.epollc
    v[ieepolhc] = ds.epolhc
    v[ieepolht] = ds.epolht
    v[ieepollt] = ds.epollt
    # ---- combustor efficiency ----
    v[ieetab]   = ds.etab
    # ---- fan map constants ----
    v[iepifK]   = ds.pifK
    v[ieepfK]   = ds.epfK
    # ---- spool speeds ----
    v[ieNf]     = eng.Nf
    v[ieN1]     = eng.N1
    v[ieN2]     = eng.N2
    v[ieNbf]    = eng.Nbf
    v[ieNblc]   = eng.Nblc
    v[ieNbhc]   = eng.Nbhc
    # ---- compressor map operating points ----
    v[iembf]    = eng.mbf
    v[iemblc]   = eng.mblc
    v[iembhc]   = eng.mbhc
    v[iepif]    = eng.pif
    v[iepilc]   = eng.pilc
    v[iepihc]   = eng.pihc
    # ---- design-point map scalars ----
    v[ieNbfD]   = ds.NbfD
    v[ieNblcD]  = ds.NblcD
    v[ieNbhcD]  = ds.NbhcD
    v[ieNbhtD]  = ds.NbhtD
    v[ieNbltD]  = ds.NbltD
    v[iembfD]   = ds.mbfD
    v[iemblcD]  = ds.mblcD
    v[iembhcD]  = ds.mbhcD
    v[iembhtD]  = ds.mbhtD
    v[iembltD]  = ds.mbltD
    v[iepifD]   = ds.pifD
    v[iepilcD]  = ds.pilcD
    v[iepihcD]  = ds.pihcD
    v[iepihtD]  = ds.pihtD
    v[iepiltD]  = ds.piltD
    # ---- duct design Mach numbers ----
    v[ieM2]     = ds.M2
    v[ieM25]    = ds.M25
    # ---- freestream (ambient) ----
    v[ieM0]     = eng.M0
    v[iep0]     = eng.p0
    v[iea0]     = eng.a0
    v[ierho0]   = eng.rho0
    v[iemu0]    = eng.mu0
    v[ieT0]     = eng.T0
    v[ieu0]     = eng.st0.u

    # ---- station total thermodynamic state (Tt, ht, pt, cpt, Rt) ----
    # st0
    v[ieTt0]   = eng.st0.Tt;  v[ieht0]  = eng.st0.ht;  v[iept0]  = eng.st0.pt
    v[iecpt0]  = eng.st0.cpt; v[ieRt0]  = eng.st0.Rt
    # st18
    v[ieTt18]  = eng.st18.Tt; v[ieht18] = eng.st18.ht; v[iept18] = eng.st18.pt
    v[iecpt18] = eng.st18.cpt;v[ieRt18] = eng.st18.Rt
    # st19
    v[ieTt19]  = eng.st19.Tt; v[ieht19] = eng.st19.ht; v[iept19] = eng.st19.pt
    v[iecpt19] = eng.st19.cpt;v[ieRt19] = eng.st19.Rt
    # st2
    v[ieTt2]   = eng.st2.Tt;  v[ieht2]  = eng.st2.ht;  v[iept2]  = eng.st2.pt
    v[iecpt2]  = eng.st2.cpt; v[ieRt2]  = eng.st2.Rt
    # st21
    v[ieTt21]  = eng.st21.Tt; v[ieht21] = eng.st21.ht; v[iept21] = eng.st21.pt
    v[iecpt21] = eng.st21.cpt;v[ieRt21] = eng.st21.Rt
    # st25
    v[ieTt25]  = eng.st25.Tt; v[ieht25] = eng.st25.ht; v[iept25] = eng.st25.pt
    v[iecpt25] = eng.st25.cpt;v[ieRt25] = eng.st25.Rt
    # st3
    v[ieTt3]   = eng.st3.Tt;  v[ieht3]  = eng.st3.ht;  v[iept3]  = eng.st3.pt
    v[iecpt3]  = eng.st3.cpt; v[ieRt3]  = eng.st3.Rt
    # st4
    v[ieTt4]   = eng.st4.Tt;  v[ieht4]  = eng.st4.ht;  v[iept4]  = eng.st4.pt
    v[iecpt4]  = eng.st4.cpt; v[ieRt4]  = eng.st4.Rt
    # st41
    v[ieTt41]  = eng.st41.Tt; v[ieht41] = eng.st41.ht; v[iept41] = eng.st41.pt
    v[iecpt41] = eng.st41.cpt;v[ieRt41] = eng.st41.Rt
    # st45
    v[ieTt45]  = eng.st45.Tt; v[ieht45] = eng.st45.ht; v[iept45] = eng.st45.pt
    v[iecpt45] = eng.st45.cpt;v[ieRt45] = eng.st45.Rt
    # st49
    v[ieTt49]  = eng.st49.Tt; v[ieht49] = eng.st49.ht; v[iept49] = eng.st49.pt
    v[iecpt49] = eng.st49.cpt;v[ieRt49] = eng.st49.Rt
    # st5
    v[ieTt5]   = eng.st5.Tt;  v[ieht5]  = eng.st5.ht;  v[iept5]  = eng.st5.pt
    v[iecpt5]  = eng.st5.cpt; v[ieRt5]  = eng.st5.Rt
    # st7 (no st6 total state in bare-pare)
    v[ieTt7]   = eng.st7.Tt;  v[ieht7]  = eng.st7.ht;  v[iept7]  = eng.st7.pt
    v[iecpt7]  = eng.st7.cpt; v[ieRt7]  = eng.st7.Rt
    # st9 — only Tt and pt in bare-pare
    v[ieTt9]   = eng.st9.Tt;  v[iept9]  = eng.st9.pt

    # ---- station static state (p, T, R, cp, u, A) ----
    # st2
    v[iep2]    = eng.st2.ps;  v[ieT2]   = eng.st2.Ts;  v[ieR2]   = eng.st2.Rs
    v[iecp2]   = eng.st2.cps; v[ieu2]   = eng.st2.u;   v[ieA2]   = ds.A2
    # st25
    v[iep25]   = eng.st25.ps; v[ieT25]  = eng.st25.Ts; v[ieR25]  = eng.st25.Rs
    v[iecp25]  = eng.st25.cps;v[ieu25]  = eng.st25.u;  v[ieA25]  = ds.A25
    # st5
    v[iep5]    = eng.st5.ps;  v[ieT5]   = eng.st5.Ts;  v[ieR5]   = eng.st5.Rs
    v[iecp5]   = eng.st5.cps; v[ieu5]   = eng.st5.u;   v[ieA5]   = ds.A5
    # st6
    v[iep6]    = eng.st6.ps;  v[ieT6]   = eng.st6.Ts;  v[ieR6]   = eng.st6.Rs
    v[iecp6]   = eng.st6.cps; v[ieu6]   = eng.st6.u;   v[ieA6]   = eng.st6.A
    # st7
    v[iep7]    = eng.st7.ps;  v[ieT7]   = eng.st7.Ts;  v[ieR7]   = eng.st7.Rs
    v[iecp7]   = eng.st7.cps; v[ieu7]   = eng.st7.u;   v[ieA7]   = ds.A7
    # st8
    v[iep8]    = eng.st8.ps;  v[ieT8]   = eng.st8.Ts;  v[ieR8]   = eng.st8.Rs
    v[iecp8]   = eng.st8.cps; v[ieu8]   = eng.st8.u;   v[ieA8]   = eng.st8.A
    # st9
    v[ieu9]    = eng.st9.u;   v[ieA9]   = eng.st9.A

    # ---- polytropic loss fractions ----
    v[ieepf]   = eng.epf
    v[ieeplc]  = eng.eplc
    v[ieephc]  = eng.ephc
    v[ieepht]  = eng.epht
    v[ieeplt]  = eng.eplt
    # ---- adiabatic efficiencies ----
    v[ieetaf]  = eng.etaf
    v[ieetalc] = eng.etalc
    v[ieetahc] = eng.etahc
    v[ieetaht] = eng.etaht
    v[ieetalt] = eng.etalt
    # ---- mass flows and power offtakes ----
    v[iemcore]  = eng.st2.mdot
    v[iemofft]  = eng.mofft
    v[iePofft]  = eng.Pofft
    v[iePhiinl] = eng.Phiinl
    v[ieKinl]   = eng.Kinl
    # ---- spool mechanical losses ----
    v[ieepsl]  = ds.epsl
    v[ieepsh]  = ds.epsh
    # ---- engine performance ----
    v[ieFe]    = eng.Fe
    v[ieFsp]   = eng.Fsp
    v[ieTSFC]  = eng.TSFC
    # ---- nozzle area schedule factors ----
    v[ieA5fac] = eng.A5fac
    v[ieA7fac] = eng.A7fac
    # ---- cooling design parameters ----
    v[iedTstrk] = ds.dTstrk
    v[ieStA]    = ds.StA
    v[ieMtexit] = ds.Mtexit
    v[ieM4a]    = ds.M4a
    v[ieruc]    = ds.ruc
    v[ieefilm]  = ds.efilm
    v[ietfilm]  = ds.tfilm
    v[iefc]     = ds.fc
    for icrow in 1:length(ds.epsrow)
        v[ieepsc1 + icrow - 1] = ds.epsrow[icrow]
        v[ieTmet1 + icrow - 1] = ds.Tmrow[icrow]
    end
    # iedeNOx (196), iePLH2 (198), ieyg (199), ieEINOx1/2 (200/201),
    # ieFAR (202), ieOPR (203), ieWc3 (204): no typed-state equivalent → NaN
    v[iemdotf]  = eng.mfuel   # total fuel mass flow (all engines)
    # ---- fuel / HX properties ----
    v[ieTft]   = eng.Tfuel_tank
    v[iehvap]  = eng.hvap
    # ieRadiatorepsilon (216), ieRadiatorMp (217): no typed-state equivalent → NaN
    v[ieRadiatorCoolantT] = eng.RadCoolantT
    v[ieRadiatorCoolantP] = eng.RadCoolantP
    v[ieRadiatorHeat]     = eng.Qradiator
    # ---- ducted-fan performance ----
    v[iemfuel]   = eng.mfuel
    v[iemfan]    = eng.mfan
    v[iePfan]    = eng.Pfan
    v[iePfanmax] = eng.Pfanmax
    v[ieTSEC]    = eng.TSEC
    # ---- HPT cooling model constants ----
    v[iefc0]     = ds.fc0
    v[iedehtdfc] = ds.dehtdfc
    # ---- convergence flag ----
    v[ieConvFail] = eng.ConvFail

    return v
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

    return deepcopy(ac.missions[imission].points[ip].engine)
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
`Dict{String,Float64}`.  Fields not written by `tfcalc!`
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
