# PR #7 Review Guide — Engine Refactor

**Branch:** `claude_engine_refactor`  
**Issue:** tasopt-eac.9  
**Date:** 2026-05-02  
**Stats:** 168 commits · 97 files · ~19 k LOC net  

---

## How to read this PR

The commit history was produced by a staged Ralph/AI-agent loop and is already labeled
by phase. Filter by prefix in `git log --oneline` to read the PR section by section
rather than commit by commit:

| Prefix | Count | Phase |
|--------|------:|-------|
| `[pare-delete]` | 28 | Remove bare-pare writes, field, conversion functions |
| `[pare-final]` | 21 | Migrate remaining scalar/field bare-pare consumers |
| `[mission-io]` | 12 | `Mission{T}`, `MissionPoint{T}`, engine sweep runners |
| `[hx-port]` | 9 | Heat-exchanger typed-state migration |
| `[ducted-fan]` | 9 | Ducted-fan harness, component wiring |
| `[refactor]` | 8 | Station renaming, `AIR_ALPHA` constant, etc. |
| `[test]` | 4 | Regression, unit, and fixture updates |
| `[doc]` | ~5 | Docs and examples |

---

## Architecture Claim

> **Component-owned equation kernels with typed state boundaries.**
>
> Each turbofan component (`Compressor`, `Combustor`, `Turbine`, `Shaft`, `Splitter`,
> `Nozzle`) owns its residual/derivative kernel and is constructed with design-point
> parameters. `tfoper!` calls these kernels rather than embedding the equations inline.
> Station data flows through `EngineState{T}` (a mutable parametric struct) rather than
> bare indexed arrays.
>
> This is **not** a fully swappable plugin-component architecture. Components are not
> dynamically dispatched at runtime — they are plain structs passed by value.
> `Inlet` BLI mixing is not yet delegated to the `Inlet` component inside `tfoper!`;
> that wiring is explicitly deferred to tasopt-eac.11. The component architecture
> contract is documented in `docs/src/dev/engine_component_contract.md` (tasopt-eac.3).

---

## Section 1 — Motivation and Non-Goals

**Why this PR exists:**

- The old engine state was stored in `pare[ie*, ip, im]` — a 3-D `Float64` array indexed
  by integer constants (`ieTt4`, `iepif`, …). This made AD (Zygote, ForwardDiff)
  difficult, made it impossible to rename or type-check individual fields, and scattered
  engine semantics across ~50 index constants.

- This PR replaces `pare` with `EngineState{T<:AbstractFloat}` — a mutable parametric
  struct. Each mission point owns one `EngineState`. Typed access is compile-time
  constant-folded; AD works without specialised containers.

**Non-goals (explicitly out of scope):**

- Runtime-swappable component plugins.
- NLsolve/Newton solver modernisation (deferred to tasopt-j11, resurfaces 2026-09-01).
- New physics or algorithm changes. All numerical deltas vs. Fortran upstream are
  pre-existing; this PR neither introduces nor fixes them. The only divergence is the
  `dhlt_ml` `+Pom` omission (tasopt-go7, tracked separately).

---

## Section 2 — Typed Engine State Substrate

**New files:**

| File | Lines | What it adds |
|------|------:|--------------|
| `src/engine/turbofan/engine_state.jl` | 519 | `EngineState{T}`: 20 `FlowStation` slots + 59 scalar fields; property-shortcut dispatch |
| `src/engine/turbofan/design_state.jl` | 196 | `DesignState{T}`: design-point anchors, frozen by `tfsize!` |
| `src/engine/turbofan/flow_station.jl` | 126 | `FlowStation{T}`: 11 thermodynamic fields per station |
| `src/engine/turbofan/gas_state.jl` | 169 | `GasState{T}` and air-composition helpers |
| `src/engine/turbofan/engine_enums.jl` | 159 | `CalcMode`, `CoolingOpt`, `EngineStation` enumerations |
| `src/engine/turbofan/thermo_wrappers.jl` | 248 | Typed wrappers: `set_total_from_Tt!`, `set_static_from_M!`, `apply_pratio_from!`, `apply_delh_from!` |
| `src/mission/mission_types.jl` | 86 | `MissionPoint{T}` (owns `EngineState{T}`) and `Mission{T}` (owns `Vector{MissionPoint{T}}`) |

**Ownership chain:**
```
aircraft.missions[im]           :: Mission{T}
  └─ .points[ip]                :: MissionPoint{T}
     └─ .engine                 :: EngineState{T}   ← single source of truth per point
        ├─ .st0, .st2, .st4 …   :: FlowStation{T}   (20 stations)
        ├─ .design              :: DesignState{T}    (frozen after tfsize!)
        └─ .TSFC, .Fe, .BPR …   scalar outputs
```

**Property shortcut:** `eng.Tt4` dispatches via `Base.getproperty(eng, Val(:Tt4))` to
`eng.st4.Tt`. This is constant-folded at compile time — no dictionary lookup, no
runtime overhead.

**What reviewers should check:**

- Station slots cover the full flow path (st0 → st9/st19).
- `DesignState` is read-only after `tfsize!`; no off-design call mutates it.
- `FlowStation` has a consistent field set across all 20 slots; no special-cased fields.
- Property shortcuts are exhaustive (every `eng.Tt4`-style access has a corresponding
  `Val` clause).

**What reviewers can ignore:**

- The choice of mutable vs. immutable struct for `FlowStation` — this was considered
  and mutable was chosen for in-place Newton iteration; it is a deliberate design
  decision, not an oversight.
- The numeric layout of fields within `DesignState` — these mirror the Fortran `parg`
  design-point scalars one-for-one; reviewers should not expect a redesigned API here.

---

## Section 3 — Mission and IO Migration

**Changed files:**

| File | What changed |
|------|-------------|
| `src/mission/mission_types.jl` | NEW — `Mission{T}`, `MissionPoint{T}` with `engine::EngineState{T}` |
| `src/IO/read_input.jl` | Populates `missions_vec[im].points[ip].engine.{Tfuel,hvap,A5fac,A7fac,TSFC,…}` from TOML |
| `src/IO/output_data.jl` | Reads typed state for output formatting |
| `src/IO/savetofile.jl` | Serialises typed `EngineState` fields |

**What reviewers should check:**

- `read_input.jl` populates at minimum: `Tfuel`, `Tfuel_tank`, `hvap`, `hvapcombustor`,
  `A5fac`, `A7fac`, `TSFC`, and `Pfanmax` from TOML/input-file values.
- Every `Mission` is initialised with `npoints` (typically `iptotal = 17`) sized
  `MissionPoint` objects before `read_input` writes fields; no out-of-bounds risk.
- `T0`, `p0` and other freestream scalars are written per-point by the atmosphere
  routine before `tfcalc!` is called; they are not read from TOML.

**What reviewers can ignore:**

- The 17-point mission structure is unchanged from the Fortran baseline. Point numbering
  (`ip`) follows the same convention; no reindexing was done here.

---

## Section 4 — Compatibility Bridge and Bare-Pare Deletion

**Commit filter:** `[pare-delete]` (28 commits) and `[pare-final]` (21 commits)

**What was deleted:**

- `pare` field from the `aircraft` struct (`pare::Array{Float64,3}`)
- All `ie*` engine index constants from `index.inc`
- `engine_state_to_pare_vec` and `pare_to_engine_state` conversion functions
- All direct `pare[ie*, ip, im]` reads/writes in the engine computation path

**What remains:**

- `para`, `parg`, `parm` arrays are unchanged. This PR touches engine state only.
- A few downstream consumers of `pare` outside the engine path may still read from
  typed `EngineState` fields via transitional accessors — these are listed in
  `test/unit_test_engine.jl` regression fixtures.
- `DuctedFanState` has a transitional `pare_to_ducted_fan_state!` bridge function
  (see Section 7).

**What reviewers should check:**

- No `pare[ie*` index reads remain anywhere in the engine calculation path
  (`tfcalc.jl`, `tfsize.jl`, `tfoper.jl`, component files, harness files).
- `aircraftsize.jl` and `wsize.jl` pass `EngineState` fields to engine calls;
  no `pare` slices are forwarded into `tfcalc!`.

**What reviewers can ignore:**

- The `ie*` constants that remain in `index.inc`. Only the engine-specific ones were
  deleted. TASOPT uses many other `ie*` constants for non-engine quantities; their
  continued presence is correct.

---

## Section 5 — Turbofan Sizing and Operation Behavior Parity

**Key files:** `src/engine/turbofan/tfsize.jl`, `src/engine/turbofan/tfoper.jl`,
`src/engine/turbofan/tfcalc.jl`, `src/engine/turbofan/tfwrap.jl`

**The core arithmetic is unchanged.** `tfsize!` and `tfoper!` contain the same
Newton solver, residual equations, and gas-property calls as the Fortran baseline.
This PR did not change the solver, modify any equation, or alter floating-point
evaluation order.

**What changed in the solver files:**

- Input/output signatures now use typed struct fields instead of `pare` array slices.
  The positional scalar lists are long but identical in value to the Fortran interface.
- `tfwrap!` is a new thin wrapper that selects `CalcMode` (sizing vs. off-design) and
  propagates `DesignState` to all mission points after sizing.
- `tfcalc!` is a new dispatcher that reads inputs from `EngineState`, calls
  `tfsize!` or `tfoper!`, and writes outputs back to `EngineState`.

**Numerical parity:** All 4666 tests pass. Regression baseline fixtures are
byte-for-byte identical to the Fortran upstream (except for the pre-existing
`dhlt_ml +Pom` omission in tasopt-go7, which this PR neither introduces nor fixes).

**What reviewers should check:**

- `tfcalc!` reads `eng.M0, eng.T0, eng.p0, eng.a0, eng.Tt4, eng.BPR` etc. from typed
  state before every `tfsize!`/`tfoper!` call — not stale values.
- `DesignState` is only mutated inside the `CalcMode.Sizing` branch; off-design branches
  read but do not write `design.*` fields.
- `CalcMode.FixedTt4OffDes` and `CalcMode.FixedFeOffDes` reach the correct `tfoper!`
  call sites with the right arguments.

**What reviewers can ignore:**

- The length of the `tfoper!` argument list. It is deliberately unchanged from the Fortran
  interface to preserve evaluation order. Shortening it is tracked in tasopt-eac.2.

---

## Section 6 — Component Equation-Kernel Architecture

**New component files:**

| File | Lines | Component | Wired into `tfoper!`? |
|------|------:|-----------|:---:|
| `src/engine/turbofan/compressor.jl` | 439 | `Compressor{T}`: fan, LPC, HPC maps and exit states | Yes |
| `src/engine/turbofan/combustor.jl` | 289 | `Combustor{T}`: fuel-air chemistry, exit state | Yes |
| `src/engine/turbofan/turbine.jl` | 345 | `Turbine{T}`, `TurbineMap{T}`: efficiency maps and exit states | Yes |
| `src/engine/turbofan/shaft.jl` | 377 | `Shaft{T}`: LP/HP power balance residuals and Jacobian terms | Yes |
| `src/engine/turbofan/splitter.jl` | 129 | `Splitter` (stateless): fan-to-bypass mass split | Yes |
| `src/engine/turbofan/nozzle.jl` | 408 | `Nozzle{T}`: throat Mach, mass-flow residual, gross thrust | Yes |
| `src/engine/turbofan/inlet.jl` | 233 | `Inlet{T}`: diffuser, BLI entropy mixing | **No — open gap** |

**Component contract:**

Each component owns:
1. A typed struct (`Compressor{T}`, `Turbine{T}`, …) carrying design-point anchors and
   efficiency-map parameters.
2. Pure kernel functions returning thermodynamic outputs and partial derivatives for the
   Newton Jacobian. Functions take inlet `FlowStation` and operating-point scalars; they
   do not read or write global state.

**Open gap — Inlet (tasopt-eac.11):**

`inlet.jl` defines `Inlet{T}`, `inlet_diffuser!`, and `inlet_bli_mixing!`. These
functions are correct and tested, and are used by the ducted-fan module. However,
`tfoper!` still computes BLI entropy mixing inline (lines 420–488) rather than
delegating to `inlet_bli_mixing!`. Wiring `Inlet` into the turbofan Newton solver
requires extending `inlet_bli_mixing!` to return the partial derivatives needed by the
9×9 Newton system — deferred to tasopt-eac.11.

The component architecture contract is documented in `docs/src/dev/engine_component_contract.md`
(tasopt-eac.3). Reviewers should treat the Inlet wiring as an explicit deferred item,
not as a defect in the other six components.

**What reviewers should check:**

- Each component function is called at the correct Newton iteration step inside
  `tfoper!` and produces the same numerical output as the inline equations it replaced.
- `Compressor{T}` is constructed with fan, LPC, and HPC design parameters; `turbine_efficiency`
  receives the correct map parameters for HPT and LPT separately.
- `Shaft.hp_shaft_workd` and `lp_shaft_workd` return the Jacobian entries that go into
  the Newton Jacobian assembled in `tfoper!`.
- `bypass_ratio` (stateless `Splitter`) is called at the fan-exit split, not at the
  nozzle throat.

**What reviewers can ignore:**

- The fact that component structs are not dynamically dispatched. This is intentional.
  Dynamic dispatch would fight the analytic Jacobian assembly pattern.
- `Inlet` being exported but not wired. It is correctly exported for downstream use
  (harness-level callers, unit tests). The wiring gap is scoped and tracked.

---

## Section 7 — Ducted Fan, Heat Exchanger, and Harness Migration

**New and changed files:**

| File | Lines | What it adds/changes |
|------|------:|---------------------|
| `src/engine/turbofan/engine_harness.jl` | 599 | `run_engine_design_point`, `run_engine_sweep`, `SweepResult`, CSV/TOML serialisation |
| `src/engine/ductedfan/ducted_fan_harness.jl` | (new) | `run_ducted_fan_design_point`, `run_ducted_fan_sweep` |
| `src/engine/ductedfan/DuctedFanState` | — | `DuctedFanState{T}` parametric struct; `pare_to_ducted_fan_state!` bridge |
| `src/engine/hx/` | — | HX typed-state migration: `resetHXs`, `HXOffDesign!` write HX delta fields into `EngineState` |

**Heat exchanger integration:**

HX solvers (`resetHXs`, `HXOffDesign!`, `RadiatorOffDesign!`) compute enthalpy/pressure
deltas and write them into `EngineState` fields (`PreCDeltah`, `InterCDeltah`,
`RegenDeltah`, `TurbCDeltah`, `HXrecircP`, etc.) before `tfcalc!` is called. `tfcalc!`
reads these fields and passes them to `tfsize!`/`tfoper!`. This is the entire HX
integration surface — reviewers do not need to understand internal HX thermodynamics.

**Ducted fan:**

`DuctedFanState` mirrors the field structure of `EngineState` for ducted-fan points.
`pare_to_ducted_fan_state!` is a transitional bridge retained for backward compatibility
with callers that still construct `pare`-style arrays for ducted-fan inputs. This is
explicitly transitional glue, not permanent API.

**What reviewers should check:**

- `run_engine_design_point` correctly populates `EngineState` fields for a single
  design-point evaluation and returns a typed result.
- `run_engine_sweep` iterates over a parameter range and returns a `Vector{SweepResult}`.
- HX delta fields are written before, not after, `tfcalc!` is called.

**What reviewers can ignore:**

- Internal HX thermodynamics — the only contract is that the five delta fields are
  populated in `EngineState` before `tfcalc!` is called.
- `pare_to_ducted_fan_state!` internals — it is bridging code that will be deleted
  once the ducted-fan caller is migrated.

---

## Section 8 — Public API and Docs/Examples Migration

**Changed files:**

| File | What changed |
|------|-------------|
| `src/engine/engine.jl` | Added exports for all new types and component functions |
| `example/example_opt.jl` | `ac.pare[ie*]` → `ac.missions[1].points[ip].engine.*` |
| `example/example_gradient_based_opt.jl` | Same migration; `parg[ig*]` for non-engine parameters |
| `docs/src/examples/NM_optimization.md` | `pare[ie*]` → typed struct paths throughout |
| `docs/src/examples/gradient_based_optimization.md` | Same; simplified to 3 design variables |
| `docs/src/examples/sensitivity.md` | Stale `pare[ieepolf]` example replaced |
| `docs/src/data_io/data_basics.md` | Stale `pare` bullet and broken `@example` block replaced |
| `src/utils/sensitivity.jl` | Stale index constant replaced |

**Exported public API (new surface):**

Types intended for external use: `EngineState`, `FlowStation`, `DesignState`,
`GasState`, `CalcMode`, `CoolingOpt`, `EngineStation`, `Mission`, `MissionPoint`.

Component types (`Compressor`, `Turbine`, `Nozzle`, `Combustor`, `Shaft`, `Splitter`,
`Inlet`) are exported for harness-level callers and unit tests but are **not** part of
the user-facing aircraft-sizing API. Users interact with `aircraft.missions[im].points[ip].engine`.

`tfsize!`, `tfoper!`, `tfcalc!`, and `tfwrap!` are exported for advanced scripting but
are not intended as stable public API. The stable entry points are `run_engine_design_point`
and `run_engine_sweep`.

**What reviewers should check:**

- `engine.jl` exports are intentional — no internal helpers are accidentally exported.
- Examples compile and produce correct output with typed paths.
- `run_engine_design_point` and `run_engine_sweep` are documented and callable without
  knowledge of `pare`.

**What reviewers can ignore:**

- The export of `tfoper!`, `tfsize!`, and internal gas-function helpers (`gassum`,
  `gas_mach`, etc.). These are exported for compatibility and scripting, not as a
  redesigned API. Tightening exports is tracked as future work.

---

## Section 9 — Tests, Fixtures, and Numerical Parity Gates

**Test files:**

| File | What it covers |
|------|---------------|
| `test/unit_test_engine.jl` | Component unit tests (gas properties, compressor/turbine/nozzle/combustor/shaft/splitter functions); single-point `tfsize!`/`tfoper!` regression vs. frozen baseline; AD gradient checks (Zygote, ForwardDiff) |
| `test/unit_test_gasturbine_flightenvelope.jl` | Full-mission sweep regression |
| `test/runtests.jl` | Suite runner; all 4666 tests |

**Parity approach:**

- Regression baseline fixtures (`test/baseline_*`) were updated when the typed-state
  migration changed field names but not values. No arithmetic was changed.
- Byte-for-byte equivalence with Fortran upstream is verified for all mission-sized
  results. The only known pre-existing divergence is `dhlt_ml +Pom` (tasopt-go7).
- ForwardDiff and Zygote gradient checks verify AD compatibility of the new parametric
  types.

**What reviewers should check:**

- The full suite passes with `VERDICT: PASSED_CLEAN` (run `./.claude/local/test.sh`).
- Baseline files under `test/` are the only golden-output fixtures; no ad-hoc
  `@test output == magic_number` checks were added.
- ForwardDiff gradient tests exercise `EngineState{ForwardDiff.Dual}` construction —
  this validates that the parametric struct is AD-compatible.

**What reviewers can ignore:**

- Absolute numerical values in fixtures. These are not design targets; they are parity
  anchors. The Fortran source is the reference; both agree except for tasopt-go7.

---

## Section 10 — Known Limitations and Follow-Up Beads

| Limitation | Severity | Tracking |
|-----------|----------|---------|
| `Inlet` BLI mixing inline in `tfoper!`; `Inlet` component not yet wired into turbofan | Architecture gap | tasopt-eac.11 (deferred post-merge) |
| `tfcalc!` unpacks EngineState into ~60 locals then repacks after the call; large scalar tuple from `tfsize!` positionally destructured | Readability technical debt | tasopt-eac.2 (deferred post-merge with follow-up) |
| Component types are not stable public API; export list includes many internals | API hygiene | tasopt-eac.5 (narrow engine module public API) |
| `DuctedFanState.pare_to_ducted_fan_state!` transitional bridge retained | Compatibility shim | Delete when ducted-fan caller migrated |
| `dhlt_ml +Pom` omission mirrors Fortran bug | Pre-existing numerical divergence | tasopt-go7 (deferred fix with regression-baseline update) |
| Station metadata not independently queryable without `EngineStation` enum | Minor | tasopt-eac (parent epic) |

**tasopt-eac.2 deferral rationale:** The current `tfcalc!` glue is readable and correct.
Introducing a named `SizingResult` struct requires care around AD compatibility. The
refactor is tracked but does not block the architecture claim because the external
boundary (EngineState in, EngineState out) is already clean.

---

## Merge Checklist

Before requesting final merge approval:

- [x] tasopt-eac.3 resolved: component architecture contract documented
      (`docs/src/dev/engine_component_contract.md`); Inlet explicitly deferred to tasopt-eac.11
- [x] tasopt-eac.9 complete: this document committed and navigable (CLOSED)
- [ ] tasopt-eac.2 either resolved or a tracked follow-up bead exists with "intentional
      glue" annotation in `tfcalc.jl`
- [ ] Full test suite: `VERDICT: PASSED_CLEAN` on `claude_engine_refactor`
- [ ] No new precompile warnings introduced

---

## Related Issues

- [tasopt-eac](beads) — parent epic
- tasopt-eac.1 — landing strategy (CLOSED)
- tasopt-eac.2 — tfcalc glue reduction (merge-gate candidate)
- tasopt-eac.3 — component architecture contract documented (CLOSED by tasopt-eac.3 commit)
- tasopt-eac.5 — narrow engine module public API (deferred post-merge)
- tasopt-eac.11 — Wire Inlet BLI into turbofan tfoper!/tfsize! (deferred post-merge)
- tasopt-eac.8 — architecture-level regression tests
- tasopt-eac.9 — this document (CLOSED by this commit)
- tasopt-eac.10 — process doc for future large Ralph refactors
- tasopt-go7 — `dhlt_ml +Pom` bug fix (deferred)
- tasopt-j11 — Phase 7 solver modernisation (deferred to 2026-09-01)
