# Turbofan Component Architecture Contract

**Issue:** tasopt-eac.3  
**Date:** 2026-05-02

---

## Overview

PR #7 introduces typed parameter structs and equation-kernel functions for the seven
core turbofan components: `Inlet`, `Compressor`, `Combustor`, `Turbine`, `Shaft`,
`Splitter`, and `Nozzle`. This document defines what the component model is, what each
component owns, and what is intentionally deferred.

---

## What a Component Is

A TASOPT turbofan component is **not** a generic pluggable object. It is a typed
**parameter container** paired with one or more **equation kernels** that operate on
those parameters and a flow state boundary. The Newton-Raphson solver in `tfoper!` owns
the iteration; components provide the physics kernels it calls.

### What a component owns

| Item | Description |
|------|-------------|
| Design-point parameters | Pressure ratio, corrected mass flow, corrected shaft speed at design |
| Map parameters | Efficiency map shape (compressors, turbines) |
| Physical constants | Fuel type, heat release, total-pressure recovery (combustor, nozzle) |
| Mechanical efficiency | Shaft gearing, bearing/windage losses |

### What a component does NOT own

- The global Newton state vector or residual/Jacobian assembly
- Station arrays or indexing into shared thermo tables
- Iteration logic or convergence checking

---

## Two-Tier API

### Tier 1 — Newton/Jacobian kernels (suffix `d`)

These are what `tfoper!` calls. Each returns a scalar or small tuple of values **plus
hand-coded partial derivatives** for the 9×9 Newton system. They are the numerically
critical path.

| Function | Component | Returns |
|----------|-----------|---------|
| `compressor_pratd` | `Compressor` | `(Nb, epol, dh, …, Jacobian terms)` |
| `combustor_burnd` | `Combustor` | `(ff, Tt4, gas composition, …)` |
| `turbine_efficiency` | `Turbine` | `(epol, ∂epol/∂dh, ∂epol/∂mb, …)` |
| `turbine_delhd` | `Turbine` | `(p_al, T_al, h_al, s_al)` — cooling-air properties |
| `turbine_mb_residual` | `Turbine` | `(residual, ∂r/∂mb)` |
| `hp_shaft_workd` | `Shaft` | `(dhht, scaling factor, Jacobian terms)` |
| `lp_shaft_workd` | `Shaft` | `(dhlt, scaling factor, Jacobian terms)` |
| `shaft_speed_residual` | `Shaft` | `(residual, ∂r/∂Nl)` |
| `bypass_ratio` | `Splitter` | `(BPR, ∂BPR/∂pt2, ∂BPR/∂pt2ac)` |
| `nozzle_exit` | `Nozzle` | `(u_e, T_e, p_e, h_e, …)` |
| `nozzle_massflow_residual` | `Nozzle` | `(residual, ∂r/∂A7)` |

All of these are wired into `tfoper!`.

### Tier 2 — State-transition kernels

These return typed `FlowStation` objects. They are designed for forward-mode simulation,
automatic differentiation, and unit-level testing — **not** for the Newton solver, which
requires explicit Jacobian terms.

| Function | Component | Purpose |
|----------|-----------|---------|
| `compressor_efficiency` | `Compressor` | Map inversion → corrected speed + polytropic efficiency |
| `compressor_exit!` | `Compressor` | Inlet FlowStation → outlet FlowStation |
| `compressor_Nb_residual` | `Compressor` | Corrected-speed map-match residual (for external Newton) |
| `combustor_exit!` | `Combustor` | Inlet FlowStation + `Tt4` → fuel fraction + outlet state |
| `turbine_exit!` | `Turbine` | Inlet FlowStation + work → outlet FlowStation |
| `hp_shaft_work` | `Shaft` | HPT specific work (scalar, no Jacobian) |
| `lp_shaft_work` | `Shaft` | LPT specific work (scalar, no Jacobian) |
| `nozzle_gross_thrust` | `Nozzle` | Momentum + pressure thrust (scalar) |
| `inlet_diffuser!` | `Inlet` | Station 0 → station 1.8 (pressure recovery) |
| `inlet_bli_mixing!` | `Inlet` | BLI entropy injection (fan face) |

Tier-2 kernels are exported from the `engine` module and covered by unit tests in
`test/unit_test_engine.jl`.

> **Why not call Tier-2 from `tfoper!`?**  
> The Newton solver needs analytic partial derivatives that Tier-2 kernels do not return.
> Rewriting Tier-1 kernels to use `FlowStation` pass-through would either require
> dual-number AD (not yet adopted for the solver path) or would duplicate the hand-coded
> Jacobian work inside the struct boundary. Tier-2 is the clean forward-simulation API;
> Tier-1 is the numerically hardened Newton API.

---

## Component Table

| Component | Struct | Tier-1 kernel(s) in `tfoper!` | Tier-2 kernels | Notes |
|-----------|--------|-------------------------------|----------------|-------|
| `Compressor` | `Compressor{T}` | `compressor_pratd` (×3: fan, LPC, HPC) | `compressor_efficiency`, `compressor_exit!`, `compressor_Nb_residual` | Map via `Ncmap`/`ecmap`; also used by ducted-fan |
| `Combustor` | `Combustor{T}` | `combustor_burnd` | `combustor_exit!` | Fuel type + LHV embedded in struct |
| `Turbine` | `Turbine{T}` | `turbine_efficiency`, `turbine_delhd`, `turbine_mb_residual` (×2: HPT, LPT) | `turbine_exit!` | |
| `Shaft` | `Shaft{T}` | `hp_shaft_workd`, `lp_shaft_workd`, `shaft_speed_residual` | `hp_shaft_work`, `lp_shaft_work` | Gear ratio and mechanical efficiency |
| `Splitter` | `Splitter` | `bypass_ratio` (×2) | — | Stateless; no design-point parameters beyond what caller provides |
| `Nozzle` | `Nozzle{T}` | `nozzle_exit`, `nozzle_massflow_residual` (×2: fan, core) | `nozzle_gross_thrust` | Area and total-pressure recovery |
| `Inlet` | `Inlet{T}` | **None — deferred (tasopt-eac.11)** | `inlet_diffuser!`, `inlet_bli_mixing!` | Used by ducted fan; turbofan still uses inline BLI |

---

## Design-Point Sizing (`tfsize!`)

`tfsize!` performs the design-point sizing pass using **direct gas function calls**
(`gas_prat`, `gas_delh`, `gas_burn`). It does not construct component instances or call
Tier-1/Tier-2 kernels. This is intentional:

- Design-point sizing solves a sequential sweep, not a coupled Newton system. It does
  not assemble a Jacobian.
- The direct gas calls are numerically equivalent to what the component kernels wrap;
  introducing typed intermediaries would not reduce duplication without a larger
  refactor.

A future refactor of `tfsize!` to share design-condition helpers with the component
definitions (e.g., a shared `design_efficiency` call) is tracked separately.

---

## Inlet Gap

The `Inlet` struct and its Tier-2 kernels (`inlet_diffuser!`, `inlet_bli_mixing!`) are
currently used only by the ducted-fan module. The turbofan solver (`tfoper!`) still
computes BLI entropy injection inline (lines 420–488 of `tfoper.jl`) with hand-coded
partials for the Newton system.

Wiring `Inlet` into `tfoper!` requires extending `inlet_bli_mixing!` to also return the
five partial derivatives needed by the Newton residual assembly (`sbfan_mf`, `sbfan_ml`,
`sbfan_Mi`, etc.). This is deferred to **tasopt-eac.11**.

---

## Extension Guidance

To add a new turbofan component:

1. Define a parameter struct `MyComp{T<:AbstractFloat}` in a new file under
   `src/engine/turbofan/`.
2. Implement at least one Tier-1 kernel returning the value plus all partial derivatives
   needed by the Newton assembly.
3. Optionally implement Tier-2 state-transition kernels returning `FlowStation` for
   forward simulation.
4. Construct the instance inside `tfoper!` from design-point parameters and call the
   Tier-1 kernel at the appropriate Newton residual line.
5. Export from `engine.jl` and cover with invariant tests in
   `test/unit_test_engine.jl`.

Do **not** move the Newton residual assembly or Jacobian accumulation into the
component. The solver structure stays in `tfoper!`; components are physics kernels only.
