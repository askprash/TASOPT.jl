# Protocol for Large Ralph/AI-Driven Refactors

**Issue:** tasopt-eac.10  
**Date:** 2026-05-02

---

## Motivation

PR #7 (`claude_engine_refactor`) was produced by a staged Ralph/AI-agent loop: 168 commits
across 97 files, ~19 k LOC net. The commit discipline was good — commits were labeled by
phase, tests passed, regression baselines were preserved. But the lead maintainer still
received one large PR whose architectural story required a separate review guide to
navigate.

This document captures the lessons. The goal is not to slow AI-assisted refactors down, but
to ensure that the **review narrative and architectural claim keep pace with the code** so
the human maintainer can act as an effective review partner, not just a final gatekeeper.

---

## Phase 0 — Refactor Charter (before the loop starts)

Before kicking off any Ralph loop that will touch more than ~10 files or ~500 LOC, file a
**charter issue** (type: `epic`) that answers:

| Question | What to write |
|----------|---------------|
| **Architecture claim** | One paragraph: what the code will look like after, and why that is better. Name the specific structural property (e.g., "typed state boundaries", "component-owned equation kernels"). |
| **Non-goals** | At least two. What this refactor will NOT do. Prevents scope creep mid-loop. |
| **Behavior invariants** | Which regression baselines must remain byte-identical, and which are allowed to change (with an upstream cross-check requirement per CLAUDE.md). |
| **Deletion boundary** | At what point is pare/legacy code removed? Before or after the new code is tested end-to-end? |
| **Exit criterion** | The specific sentence that will close this epic ("All 6 typed components are wired into `tfoper!` and the inline duplicate is removed"). |

Tag the charter `architecture` and link it from any `docs/src/dev/` doc it affects.

---

## Phase 1 — Invariant Tests Before Migration

Before modifying any production file, write tests that assert the **structural invariants
the refactor is supposed to achieve**. These are not golden-value snapshot tests; they are
property tests that will fail if the refactor drifts from its claim.

Examples:
- "Every exported engine component struct has a parametric type `T<:AbstractFloat`."
- "The set of fields in `EngineState` matches `_STATION_DUMP_ORDER` exactly."
- "No symbol in the engine module's public export list starts with `_`."

These tests form a **behavior gate** — they must pass before the epic can be closed, and
they protect the architectural claim against future regression.

Commit these tests on their own before the loop proceeds. If the loop runs in parallel,
the test commit is the loop's first task.

---

## Phase 2 — Living Review Guide

As soon as the first commit lands, create `docs/src/dev/<name>_review_guide.md`. Update it
in **the same Ralph batch that produces each major phase commit**. The guide must always be
current enough that a reviewer could pick it up and understand what just changed.

Minimum sections:

```
## Architecture Claim
## Phase Map  (label prefix → what was being done)
## Open Gaps  (what is not yet done; link bd issue for each)
## Merge Checklist
```

The review guide is not a retrospective artifact — it is written in parallel with the code
so that the maintainer can review incrementally rather than all at once at the end.

If the guide is more than one week stale relative to the latest commit, treat that as a
process failure and file a bd issue to update it before continuing.

---

## Phase 3 — Review Checkpoints Before Irreversible Steps

Some refactor steps are easy to reverse; others are not. Before any **irreversible
deletion phase** (removing a legacy struct, dropping a module, deleting compatibility
glue), pause and satisfy all three:

1. **Tests pass end-to-end** on the post-deletion state.
2. **Review guide is current** — open gaps are documented, not silently hidden.
3. **Human sign-off** if the maintainer is not also the agent operator. File a
   `review-checkpoint` issue; close it only after the maintainer has read the current guide
   and explicitly OKs the deletion.

"Irreversible" means: a reviewer cannot reconstruct the deleted code from git history
without significant effort (i.e., the deletion is non-trivially entangled with other
commits).

---

## Phase 4 — Compatibility Glue Tagging

Any code that exists **only** to bridge old callers to new internals is compatibility glue.
Examples: a wrapper that reads from `EngineState` and writes back to a bare array for a
single caller; a type alias that re-exports a renamed struct.

Rules for glue:

1. **Mark it.** Add a comment `# COMPAT(bd-id): reason` at the glue site. The bd issue
   must exist before the glue is committed — not after.
2. **Scope it.** Each glue site has exactly one bd issue. Do not create unbounded glue
   with a vague "clean up later" comment.
3. **Delete it on schedule.** The bd issue for glue removal must be filed in the same
   commit that introduces the glue, set to the same priority as the feature it enables.
   Glue that has no deletion issue is permanent by default, even if it was intended to be
   temporary.

---

## Phase 5 — Stacked PRs vs. One PR with Sections

### Use stacked PRs when

- The refactor has a **clean plateau**: a commit (or small set of commits) where tests
  pass, baselines are clean, and the architecture is coherent — even if incomplete.
- Each stack level can be reviewed in under ~2 hours by someone who knows the codebase.
- The phases are **not interleaved**: phase A commits and phase B commits do not need to
  be cherry-picked around each other.

### Use one PR with review sections when

- Phases are interleaved in commit history (as in PR #7, where `[pare-delete]` and
  `[pare-final]` commitsare entangled).
- Retro-splitting would require re-creating intermediate compatibility states that no
  longer exist cleanly. The cost of the split (4–8 hours of rebase work) exceeds the
  review burden.
- The maintainer is also the primary reviewer and prefers a single coherent narrative
  over multiple stacked review rounds.

In the one-PR case, the review guide substitutes for a clean commit boundary: each section
of the guide corresponds to a filter on the commit prefix (e.g., `[pare-delete]` commits,
then `[mission-io]` commits). The guide must make this explicit.

**Decision threshold:** if a PR requires more than one sitting to review (roughly > 4 hours
for a familiar reviewer), it should be stacked if the phases are separable. If they are
not separable, the guide must be good enough to allow incremental review across multiple
sittings.

---

## Phase 6 — Behavior Gates

At minimum, two behavior checkpoints must exist in any large refactor:

| Checkpoint | When | What |
|------------|------|------|
| **Pre-migration baseline** | Before the first production change | Run full suite; record test count and baseline hash. Store in progress.txt. |
| **Post-migration gate** | Before closing the epic | Full suite passes; baselines byte-identical (or delta is justified per CLAUDE.md). Invariant tests from Phase 1 all pass. |

If the refactor spans multiple Ralph sessions, the post-session note in `progress.txt`
must record the current baseline status ("PASSED_CLEAN N/N, baselines unchanged" or "N
tests, delta in X justified by tasopt-go7").

---

## Checklist for Starting a New Large Refactor

```
[ ] Charter issue (epic) filed with architecture claim, non-goals, and exit criterion
[ ] Invariant tests committed before first production change
[ ] Review guide created (docs/src/dev/<name>_review_guide.md)
[ ] Pre-migration baseline recorded in progress.txt
[ ] Deletion phase checkpoint issue filed if any legacy code will be removed
[ ] All compatibility glue sites have COMPAT(bd-id) tags and companion deletion issues
[ ] Decision recorded: stacked PRs or one PR with sections (with rationale)
```

---

## Relationship to CLAUDE.md Refactor/Logic Rules

This protocol is layered on top of, not a replacement for, the CLAUDE.md rules:

- **Refactor vs. logic change** (CLAUDE.md §Refactor vs Logic Change) still applies at
  the commit level. This protocol operates at the epic/PR level.
- The "behavior gate" checkpoints here enforce the CLAUDE.md rule that any bit-level delta
  in baselines requires its own bd issue with upstream cross-check.
- "Compatibility glue tagging" is an extension of the anti-pattern documented in
  CLAUDE.md §Examples (e83166a7): glue that is not explicitly tracked becomes a silent
  change.

---

## Acknowledgment

Some refactors are large because the technical debt is real, not because the process
failed. A well-structured agent loop that produces 168 labeled commits is better than a
human author who delivers one monolithic diff. The protocol above does not make large
refactors smaller — it keeps the review narrative honest and the human maintainer in the
loop at each irreversible step.
