"""
    generate_engine_state_baseline.jl

Standalone script to regenerate the compact engine-state regression baseline.

## Usage

    julia --project=. test/generate_engine_state_baseline.jl

## When to regenerate

Run this script when a deliberate physics or model change causes the engine
state to shift from the pinned reference values.  The updated baseline must
be committed together with the code change so that reviewers can audit the
diff side-by-side.

## What it produces

    test/fixtures/default_sized_engine_state.toml

A compact TOML file (~30 lines) containing:
- `full_state_hash`: SHA-256 of the canonicalized full engine state (all
  mission points, same fields as the old per-field loop, floats formatted
  as "%.10e", keys sorted as "ip{ip:02d}/path").
- `[cruise1]`: 8 curated endpoint values at ipcruise1 that an engineer
  can reason about when auditing a regression diff.

See test/fixtures/README.md for the full canonicalization spec.
"""

using TASOPT

@warn """
*** REGENERATING ENGINE STATE BASELINE ***
This script rewrites test/fixtures/default_sized_engine_state.toml.
Run it only when a deliberate physics or model change is intended.
Commit the updated baseline with an explanatory git message that includes
a BASELINE-REGEN: marker and an upstream cross-check note (per CLAUDE.md).
"""

ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)
path = TASOPT.reset_regression_test_engine_state(ac)
println("Wrote $path")
