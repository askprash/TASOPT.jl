# Project Context

This is **`TASOPT.jl`** — the main working directory (Julia port of TASOPT). All code changes and new work happen here.

The sibling folder **`../Tasopt2.16/`** is the original Fortran TASOPT 2.16 implementation, kept as **reference only**. Consult it when you need to understand original algorithms or numerical behavior, but **do not modify it**.

## Refactor vs Logic Change

A **refactor** commit changes structure only: extract a function, rename a variable, wrap values in a typed container, reorder declarations. The arithmetic, branch structure, and floating-point evaluation order must remain identical. A clean refactor produces **zero diff** in regression baselines.

A **logic change** alters behavior: different algorithm, different branch structure, different FP operation order, new physics, a bug fix. Any commit that changes a regression baseline — even by 1 ULP — is a logic change.

### Rules

1. **Never bundle a logic change inside a refactor commit.** If you discover a bug while refactoring, file a bd issue for the bug and fix it in a separate commit with its own regression-baseline update.
2. **Any bit-level delta in regression baselines requires its own bd issue** with an upstream cross-check explaining why the delta is correct.
3. **FP evaluation order matters.** Rewriting `a * (b + c)` as `a*b + a*c` is a logic change, not a refactor, because it changes the rounding sequence. Extracting `a * (b + c)` into a helper function that evaluates the same expression is a refactor.

### Examples

- **Anti-pattern (e83166a7):** Commit `e83166a7` (tasopt-j9l.58) extracted `lp_shaft_workd` as a Jacobian wrapper (refactor) but simultaneously corrected a missing `+Pom` term in `dhlt_ml` (bug fix). The bug fix silently changed regression baselines. These should have been two separate commits: one for the wrapper extraction, one for the Pom correction with its own bd issue (tasopt-go7).
- **Good pattern (j9l.27, shaft wrappers):** Wrapper extraction preserved upstream arithmetic exactly — byte-for-byte baseline match confirmed.
- **Good pattern (j9l.55, splitter):** Self-corrected back to byte-for-byte equivalence with upstream after initially drifting.


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
