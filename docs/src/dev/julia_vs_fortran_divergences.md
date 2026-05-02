# Julia-vs-Fortran Known Divergences

This page records every intentional or accepted divergence between the Julia port
(`TASOPT.jl`) and the upstream Fortran reference (`TASOPT 2.16`).  Each entry must
be added when the divergence is introduced — not retroactively — so that reviewers
can understand the delta in any pull request.

## How to add an entry

Copy the template below, fill in every field, and append it under the appropriate
status section.  Entries are ordered by date introduced within each section.

```markdown
### [Short title]

| Field                | Value |
|----------------------|-------|
| **Julia location**   | `src/path/to/file.jl`, function `foo`, lines N–M |
| **Fortran location** | `src/path/to/file.f`, lines N–M |
| **Nature**           | One- or two-sentence summary of what differs. |
| **Why Julia is correct / intentional** | Explanation of correctness or intent. |
| **Steady-state delta** | Quantified impact at the converged fixed point, e.g. `~1.56e-8` relative in key output, or "none measurable at converged fixed point." |
| **Owning bd issue**  | `tasopt-XXXX` |
| **Date introduced**  | YYYY-MM-DD |
```

---

## Pending / deferred fixes

Entries in this section track divergences that are *planned* but not yet landed —
typically deferred bug fixes that will be applied in a future iteration.

### LP-shaft dhlt\_ml partial missing `+Pom` in demand factor

| Field                | Value |
|----------------------|-------|
| **Julia location**   | `src/engine/turbofan/tfoper.jl`, computation of `dhlt_ml`, lines 1372–1374 |
| **Fortran location** | `src/engine/turbofan/tfoper.f` (Fortran reference copy), line ~1423 |
| **Nature**           | `dhlt_ml` computes the LP-shaft-work partial w.r.t. corrected core-fan mass flow using `demand_no_pom = ht2_5 - ht2ac + BPR*(ht13 - ht2)`, omitting the `+Pom` term that appears in every other `dhlt_*` partial.  This is a copy-paste omission present in both the Fortran and the current Julia code. |
| **Why Julia is correct / intentional** | Once `tasopt-go7` is applied, Julia will add `+Pom` to the `dhlt_ml` demand factor, matching the pattern of all other `dhlt_*` partials.  The Fortran retains the bug.  The error magnitude is `Pom * dlfac_ml`, giving ~1.56e-8 relative difference in the converged solution; the primary impact is a slightly inexact Jacobian that perturbs Newton convergence. |
| **Steady-state delta** | ~1.56e-8 relative in converged solution variables (quantified during audit tasopt-9r4 / tasopt-cv4).  Baselines will be updated when the fix lands. |
| **Owning bd issue**  | `tasopt-go7` |
| **Date introduced**  | 2026-09-01 (planned landing date; entry pre-filed 2026-05-01) |

---

## Active divergences

*No entries yet.*  The first real entry will be added when `tasopt-go7` lands.
