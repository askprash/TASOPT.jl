# test/fixtures — Regression Baselines

## Files

| File | Testset | Regen script |
|------|---------|--------------|
| `default_sized_engine_state.toml` | `regression_test_size_aircraft.jl` → "Propulsion (typed EngineState)" | `test/generate_engine_state_baseline.jl` |
| `engine_sweep_baseline.toml` | `unit_test_engine.jl` → "Engine sweep" | `test/generate_engine_baseline.jl` |
| `ducted_fan_sweep_baseline.toml` | `unit_test_ductedfan.jl` | `test/generate_ducted_fan_baseline.jl` |
| `engine_benchmark_baseline.toml` | `unit_test_engine.jl` → "Engine benchmark" | `test/generate_engine_benchmark_baseline.jl` |

---

## default_sized_engine_state.toml — canonicalization scheme

The file stores two things:

### 1. `full_state_hash` — SHA-256 over the full engine state

Covers all `nip` mission points (17 for the default aircraft) and the same
field set as the old per-field assertion loop:

- **EngineState scalars** — `M0`, `T0`, `p0`, `Fe`, `TSFC`, `BPR`, … (53 fields)
- **DesignState scalars** — `pi_fan_des`, `mb_fan_des`, `Nb_fan_des`, … (47 fields)
- **DesignState vectors** — `epsrow[1:4]`, `Tmrow[1:4]` (8 elements)
- **FlowStation fields** at all 20 stations — `Tt`, `ht`, `pt`, `cpt`, `Rt`, `Ts`, `ps`, `cps`, `Rs`, `u`, `A`, `mdot` (240 values per point)

**Step-by-step**:

1. For each mission point `ip = 1 .. nip`, for each value in the above set,
   form a key `"ip{ip:02d}/{path}"` (e.g. `"ip10/Fe"`, `"ip10/st4/Tt"`).
2. Sort all (key, value) pairs lexicographically by key.
3. Format each float as `"%.10e"` (11 significant figures). This is well above
   cross-platform libm ULP noise (~10⁻¹⁵) but tight enough to catch any real
   physics regression (changes of order 10⁻¹⁰ or larger).
4. Write one `"key=value\n"` line per pair to an in-memory buffer.
5. Compute SHA-256 of the resulting byte string and encode as lowercase hex.

The canonicalization function lives in two places that **must be kept in sync**:
- `_engine_state_canonical_hash` in `test/regression_test_size_aircraft.jl` (used by the test)
- `_engine_state_canonical_hash_src` in `src/IO/save_model.jl` (used by `reset_regression_test_engine_state`)

### 2. `[cruise1]` — curated endpoint spot-checks

Eight physically meaningful values at `ipcruise1` (ip = 10 for the default
mission):

| Key | Field | Units | Description |
|-----|-------|-------|-------------|
| `Fe` | `eng.Fe` | N | Net thrust |
| `TSFC` | `eng.TSFC` | kg/(N·s) | Thrust-specific fuel consumption |
| `BPR` | `eng.BPR` | — | Bypass ratio |
| `Tt4` | `eng.st4.Tt` | K | Turbine-inlet total temperature |
| `pt3` | `eng.st3.pt` | Pa | Compressor-exit total pressure |
| `Nbf` | `eng.Nbf` | — | Fan corrected speed |
| `Nblc` | `eng.Nblc` | — | LPC corrected speed |
| `Nbhc` | `eng.Nbhc` | — | HPC corrected speed |

These are checked with `rtol=1e-10`, the same tolerance used in the old
per-field loop.

---

## When to regenerate

Regenerate **only** when a deliberate physics or model change causes the
engine state to shift. Per `CLAUDE.md`:

- Any change to a regression baseline is a **logic change**, not a refactor.
- The commit must include a `BASELINE-REGEN:` marker and an upstream
  cross-check note explaining why the delta is correct.

### Steps

```bash
julia --project=. test/generate_engine_state_baseline.jl
# Review the diff — hash changed + which curated values changed
git diff test/fixtures/default_sized_engine_state.toml
# Commit with BASELINE-REGEN: marker
git add test/fixtures/default_sized_engine_state.toml
git commit -m "[fix|refactor]: <description>

BASELINE-REGEN: <explain why the delta is correct, cross-check against Fortran>"
```

### Verifying a synthetic perturbation

To confirm the test catches real regressions, edit one numeric value in the
TOML (e.g. change `Fe` by 1 ULP) and run the test suite — you should see
the hash assertion fail with a clear mismatch message.
