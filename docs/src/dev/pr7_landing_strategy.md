# PR #7 Landing Strategy

**Issue:** tasopt-eac.1  
**Date:** 2026-05-02  
**Decision:** Keep as **one PR** with a hardened review guide. Do not retro-split.

---

## Context

PR #7 (`claude_engine_refactor`) contains 168 commits across 97 files (~19 k LOC net).
The commit history was produced by a staged Ralph/AI-agent loop and is already labeled
by phase:

| Prefix | Count | Content |
|--------|------:|---------|
| `[pare-delete]` | 28 | Remove bare-pare writes, field, conversion functions |
| `[pare-final]` | 21 | Migrate remaining scalar/field bare-pare consumers |
| `[mission-io]` | 12 | `Mission{T}`, `MissionPoint{T}`, engine sweep runners |
| `[hx-port]` | 9 | Heat-exchanger typed-state migration |
| `[ducted-fan]` | 9 | Ducted-fan harness, component wiring |
| `[refactor]` | 8 | Station renaming, `AIR_ALPHA` constant, etc. |
| `[test]` | 4 | Regression, unit, and fixture updates |
| `[doc]` | ~5 | Docs and examples |

New files: `engine_state.jl`, `gas_state.jl`, `design_state.jl`, `flow_station.jl`,
`thermo_wrappers.jl`, `engine_harness.jl`, `ducted_fan_harness.jl`, `mission_types.jl`,
plus seven component files (`inlet.jl`, `compressor.jl`, `combustor.jl`, `turbine.jl`,
`nozzle.jl`, `shaft.jl`, `splitter.jl`).

---

## Decision: One PR — No Retro-Split

### Rationale

1. **Commit history is already staged.** A reviewer can filter by prefix to read the PR
   in phases without needing branch splits.

2. **Retro-splitting is prohibitively expensive.** The `[pare-delete]` and `[pare-final]`
   phases are interleaved with component wiring. Cherry-picking 168 commits across three
   stacked PRs would require re-creating intermediate compatibility states that no longer
   exist cleanly.

3. **Tests pass end-to-end.** All 4666 tests pass; regression baselines are byte-for-byte
   identical to Fortran upstream (where applicable). A split would need to pass tests at
   each stage, which forces intermediate compatibility that adds complexity without value.

4. **The maintainer is also the primary reviewer.** A well-structured review guide
   (tasopt-eac.9) is a more efficient substitute for branch splitting.

5. **The architecture claim is accurate enough.** Six of seven component types
   (`Compressor`, `Combustor`, `Turbine`, `Shaft`, `Splitter`, `Nozzle`) are wired into
   `tfoper!` via equation-kernel methods. `Inlet` BLI mixing remains inline in `tfoper!`
   — this is the one gap (tracked in tasopt-eac.3).

---

## Gates Before Merge

These must be resolved **before** this PR can be merged. They are tracked as open child
issues of tasopt-eac:

| Issue | Title | Status |
|-------|-------|--------|
| **tasopt-eac.3** | Make turbofan component extraction coherent | OPEN — must resolve or explicitly narrow Inlet claim |
| **tasopt-eac.9** | Write PR #7 review guide and architecture claim | OPEN — must complete before requesting final review |
| **tasopt-eac.2** | Reduce tfcalc typed-state projection glue | OPEN — resolve or file a tracked follow-up if deferring |

`tasopt-eac.2` (tfcalc glue) may be deferred post-merge if a follow-up bead is filed
and the current state is documented as intentional glue. `tasopt-eac.3` and `tasopt-eac.9`
are hard gates: the PR cannot claim a modular component architecture without resolving the
Inlet gap, and it cannot request final review without a review guide.

---

## Fallback: Minimal Two-Way Split

If, after the review guide is written, the PR is still judged unreviewable in one pass,
the **minimum viable split** is at one boundary:

- **PR A:** Typed-state substrate + pare deletion
  (`EngineState`/`GasState`/`DesignState` introduction through final `[pare-delete]`
  commits; no component types; all regression tests pass)
- **PR B:** Component extraction + ducted-fan/HX/mission-IO + docs/API/tests
  (built on top of PR A; claims component equation-kernel architecture)

This split is feasible because the pare-deletion phase forms a coherent plateau (typed
state present, bare pare gone, tests pass). Component wiring builds on top of that
plateau without backtracking.

Only pursue this fallback if the review guide is insufficient. The cost of re-creating
the split point cleanly is estimated at 4–8 hours of rebase work.

---

## Architecture Claim

The architecture claim for PR #7 is:

> **Component-owned equation kernels with typed state boundaries.**
>
> Each turbofan component (`Compressor`, `Combustor`, `Turbine`, `Shaft`, `Splitter`,
> `Nozzle`) owns its residual/derivative kernel and is constructed with design-point
> parameters. `tfoper!` calls these kernels rather than embedding the equations inline.
> `Inlet` BLI mixing is not yet delegated to the `Inlet` component in `tfoper!` — this
> is an open gap tracked in tasopt-eac.3.
>
> This is **not** a fully swappable plugin-component architecture. Components are not
> dynamically dispatched at runtime. The typed state boundary means station data flows
> through `EngineState` rather than bare indexed arrays.

---

## Related Issues

- [tasopt-eac](https://github.com/MIT-LAE/TASOPT.jl) — parent epic
- tasopt-eac.1 — this decision (CLOSED on completion)
- tasopt-eac.2 — tfcalc glue reduction
- tasopt-eac.3 — component architecture coherence (Inlet gap)
- tasopt-eac.9 — PR review guide
- tasopt-eac.10 — process doc for future large Ralph refactors
