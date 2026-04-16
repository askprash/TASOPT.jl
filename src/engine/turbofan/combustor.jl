"""
Combustor component type for a TASOPT turbofan.

Models the main burner from the HPC exit (station 3) to the turbine inlet
(station 4), computing the fuel-to-air ratio and the exit gas composition
for a given outlet total temperature.
"""

# ---------------------------------------------------------------------------
# Combustor struct
# ---------------------------------------------------------------------------

"""
    Combustor{T<:AbstractFloat}

Parametric combustor component holding burner pressure-recovery, combustion
efficiency, and fuel-stream properties.

Parametric in the numeric type `T` so that forward-mode AD and other dual-
number types propagate without requiring specialised containers.  Note that
`ifuel` is an `Int` index (not `T`) because it selects a fuel species and
carries no numerical sensitivity.

## Fields

| Field   | Unit  | Description                                                     |
|:--------|:------|:----------------------------------------------------------------|
| `pib`   | —     | Burner total-pressure ratio `pt4 / pt3` (< 1, typ. 0.94–0.98)  |
| `etab`  | —     | Combustion efficiency — fraction of fuel that is fully burned;  |
|         |       | `1 − etab` remains as unburnt fuel in the product stream        |
| `Ttf`   | K     | Fuel total temperature at the burner inlet                      |
| `ifuel` | —     | Integer index identifying the fuel species; passed to `gasfuel` |
|         |       | and `gasfun` for combustion-product chemistry                   |
| `hvap`  | J/kg  | Fuel latent heat of vaporisation; subtracted from fuel enthalpy |
|         |       | in `gas_burnd` before computing the energy balance              |

## Physics summary

Given the compressor-exit total state (station 3) and a target burner-exit
total temperature `Tb` (= `Tt4`), `combustor_exit!` solves the enthalpy
balance

    ṁ_air · Δh_air  =  ṁ_fuel · (h_fuel − h_products)

via `gas_burnd`, yielding the fuel-to-burner-inlet-flow ratio `ffb` and the
mixed exit composition `lambda`.  The exit total pressure is then simply

    pt4 = pib · pt3

and the exit thermodynamic state follows from `gassumd(lambda, nair, Tt4)`.

## Station convention

| Station | Description                             |
|:--------|:----------------------------------------|
| 3       | HPC exit / burner inlet                 |
| 4       | Burner exit (before cooling-air mixing) |
"""
mutable struct Combustor{T<:AbstractFloat}
    pib  ::T    # burner total-pressure ratio  pt4/pt3      [—]
    etab ::T    # combustion efficiency                      [—]
    Ttf  ::T    # fuel inlet total temperature               [K]
    ifuel::Int  # fuel species index (integer, not T)
    hvap ::T    # fuel enthalpy of vaporisation              [J/kg]
end

# ---------------------------------------------------------------------------
# combustor_exit!  —  energy balance → fuel fraction + exit state
# ---------------------------------------------------------------------------

"""
    combustor_exit!(st4, st3, burner, alpha, Tb)
        -> (ffb, lambda, ffb_Tt3, ffb_Ttf, ffb_Tb, lam_Tt3, lam_Ttf, lam_Tb)

Compute the burner-exit `FlowStation` and fuel-to-air ratio for the combustor
`burner`, given the inlet state `st3` and the target exit total temperature
`Tb`.

## Arguments

| Argument | Unit  | Description                                                   |
|:---------|:------|:--------------------------------------------------------------|
| `st4`    | —     | Outlet `FlowStation` (mutated in place)                       |
| `st3`    | —     | Inlet `FlowStation` (read-only); must have valid `pt`, `Tt`,  |
|          |       | `ht`, `st`, `cpt`, `Rt`, and `alpha` fields                   |
| `burner` | —     | `Combustor` component                                         |
| `alpha`  | —     | Air mass-fraction vector (length `n = 6`); typically the      |
|          |       | standard atmospheric composition used throughout `tfoper!`    |
| `Tb`     | K     | Target burner-exit total temperature (= `Tt4`)                |

## Returns

An 8-tuple:

    (ffb, lambda, ffb_Tt3, ffb_Ttf, ffb_Tb, lam_Tt3, lam_Ttf, lam_Tb)

| Return      | Description                                              |
|:------------|:---------------------------------------------------------|
| `ffb`       | Fuel-to-burner-inlet-flow mass ratio `ṁ_fuel / ṁ_burner`|
| `lambda`    | Exit gas composition vector (length `n`)                 |
| `ffb_Tt3`   | ∂ffb/∂Tt3                                                |
| `ffb_Ttf`   | ∂ffb/∂Ttf                                                |
| `ffb_Tb`    | ∂ffb/∂Tb                                                 |
| `lam_Tt3`   | ∂lambda/∂Tt3 (length-n vector)                           |
| `lam_Ttf`   | ∂lambda/∂Ttf (length-n vector)                           |
| `lam_Tb`    | ∂lambda/∂Tb  (length-n vector)                           |

The derivative arrays are needed by the Newton driver (tasopt-j9l.31) for
Jacobian assembly.

## Physics

1. **Combustion-product chemistry** — `gasfuel(ifuel, n)` returns the
   per-species mass-fraction changes `gamma` due to complete combustion.
   Partial efficiency is applied: `gamma[i] *= etab` for air species, and
   `gamma[n] = 1 − etab` accounts for unburnt fuel in the product stream.

2. **Energy balance** — `gas_burnd` solves

       ṁ_air·(h(Tb) − h(Tt3)) = ffb·(h_fuel(Ttf) − h_products(Tb))

   returning `ffb` and the mixed exit composition `lambda`.

3. **Pressure recovery** — `pt4 = pib · pt3` (no work, simple loss).

4. **Exit thermodynamic state** — `gassumd(lambda, nair, Tt4)` populates
   `st4.st`, `st4.ht`, `st4.cpt`, `st4.Rt`.

## Notes

- `beta = [0,0,0,0,0,1]` is the pure-fuel stream composition and is
  constructed internally — it does not depend on the operating point.
- Static-state fields (`Ts`, `ps`, `u`, …) are **not** set by this function;
  those depend on the downstream Mach number and are set by
  `set_static_from_M!` when the flow is matched.
- `n = 6` and `nair = 5` are the TASOPT standard gas-species counts.
"""
function combustor_exit!(
    st4    ::FlowStation{T},
    st3    ::FlowStation{T},
    burner ::Combustor{T},
    alpha  ::AbstractVector,
    Tb     ::T,
) where {T<:AbstractFloat}

    # gas_burnd convention: n = nair = number of air species (5 for TASOPT).
    # The fuel species is handled internally by gas_burnd via `ifuel` and
    # `gasfun`, NOT via the alpha/beta/gamma vectors at index n+1.
    nair = 5

    # Pure-fuel stream beta: all air-species slots are zero; gas_burnd adds
    # the fuel contribution itself via gasfun(ifuel, Ttf).
    beta = zeros(T, nair)

    # --- combustion-product chemistry -----------------------------------
    # gasfuel returns nair-species delta fractions (changes per unit fuel burned).
    # gasfuel always operates on the 6-species layout internally; we request
    # 6 elements and use only the first nair=5 (air species).
    gamma_raw = gasfuel(burner.ifuel, nair + 1)

    # Promote to T so ForwardDiff dual numbers propagate: T(constant) gives
    # a dual with value=constant and zero partial, which is correct.
    gamma = T.(gamma_raw[1:nair])

    # Apply combustor efficiency: scale the reaction change fractions.
    # Unburnt fuel (1 − etab) is tracked implicitly by the reduced gamma
    # but is not carried as a separate species in the nair-species model.
    for i in 1:nair
        gamma[i] = burner.etab * gamma[i]
    end

    # --- energy balance: fuel fraction + exit composition ---------------
    ffb, lambda,
    ffb_Tt3, ffb_Ttf, ffb_Tb,
    lam_Tt3, lam_Ttf, lam_Tb = gas_burnd(
        alpha, beta, gamma, nair, burner.ifuel,
        st3.Tt, burner.Ttf, Tb, burner.hvap,
    )

    # --- outlet total pressure ------------------------------------------
    st4.pt = burner.pib * st3.pt

    # --- outlet total temperature and composition -----------------------
    # lambda is nair-element; FlowStation.alpha is SVector{5} — compatible.
    st4.Tt    = Tb
    st4.alpha = lambda

    # --- outlet thermodynamic state from exit composition ---------------
    # gassumd(lambda, nair, Tb) → (s, s_T, h, h_T, cp, R) at temperature Tb
    st4.st, _, st4.ht, _, st4.cpt, st4.Rt = gassumd(lambda, nair, Tb)

    return ffb, lambda, ffb_Tt3, ffb_Ttf, ffb_Tb, lam_Tt3, lam_Ttf, lam_Tb
end
