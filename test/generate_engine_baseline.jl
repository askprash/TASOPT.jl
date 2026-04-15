"""
    generate_engine_baseline.jl

Standalone script to regenerate the engine-sweep regression baseline.

## Usage

    julia --project=. test/generate_engine_baseline.jl

## When to regenerate

Run this script when a deliberate physics or model change causes the
engine sweep outputs to shift from the pinned reference values.  The
updated baseline must be committed together with the code change so that
reviewers can audit the numerical diff side-by-side.

## What it produces

    test/fixtures/engine_sweep_baseline.toml

A TOML file with one [[points]] entry per mission point (all 16 regular
points by default).  Each entry records performance scalars and the full
station-level thermodynamic state, serialised at full Float64 precision.
"""

using TASOPT

# ------------------------------------------------------------------
# Auditable warning — the script echoes the same warning that
# regenerate_engine_baseline() emits internally, so it is visible
# in CI logs if this script is run accidentally.
# ------------------------------------------------------------------
@warn """
*** REGENERATING ENGINE SWEEP BASELINE ***
This script rewrites test/fixtures/engine_sweep_baseline.toml.
Run it only when a deliberate physics or model change is intended.
Commit the updated baseline with an explanatory git message.
"""

ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)
TASOPT.engine.regenerate_engine_baseline(ac)
