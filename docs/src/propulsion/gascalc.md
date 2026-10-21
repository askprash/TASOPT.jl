# Gas Calculations

Gas calculations used with TASOPT.jl propulsion systems.

## Property data

Species properties come from [IdealGasThermo.jl](https://github.com/MIT-LAE/IdealGasThermo.jl).
Each species is an `IdealGasThermo.FrozenGas` built once at load time from the
NASA-9 polynomial fits in IdealGasThermo's database. Evaluating ``c_p``, ``h``
and ``s`` is then a single polynomial evaluation, and ``dc_p/dT`` is exact
(the earlier table-based implementation used a finite difference).

Gases are still addressed by the integer index `igas` used throughout the
engine code:

| `igas` | Species | Data |
|-------:|---------|------|
| 1  | N₂   | NASA-9 |
| 2  | O₂   | NASA-9 |
| 3  | CO₂  | NASA-9 |
| 4  | H₂O  | NASA-9 |
| 5  | Ar   | NASA-9 |
| 11 | CH₄  | NASA-9 |
| 12 | C₂H₆ | NASA-9 |
| 13 | C₃H₈ | NASA-9 |
| 14 | C₄H₁₀ (n-butane) | NASA-9 |
| 18 | C₈H₁₈ (n-octane) | NASA-9 |
| 24 | C₁₄H₃₀ (Jet-A surrogate) | NASA-9 Jet-A(g) temperature dependence, C₁₄H₃₀ formation enthalpy and molecular weight |
| 40 | H₂   | NASA-9 above 200 K, tabulated LH₂-range data below 200 K |

Notes:

- **Datums.** Enthalpy includes the heat of formation (``h(298.15\,\mathrm{K})``
  is the formation enthalpy). Entropy ``s`` is the standard-state entropy
  function ``\int c_p/T\,dT``. Only differences of ``s`` at fixed composition
  enter the cycle calculations.
- **Temperature range.** The NASA-9 fits are valid from 200 K to 6000 K. Below
  200 K, H₂ (index 40) switches to tabulated data down to 20.4 K, shifted so
  that ``h`` and ``s`` are continuous at 200 K. This range is used for LH₂
  coolant in the heat exchangers and for cold H₂ fuel. Other species are
  extrapolated below 200 K.
- **Fuel index 24.** IdealGasThermo has no C₁₄H₃₀. Index 24 uses the mass-specific
  ``c_p(T)`` and sensible enthalpy of IdealGasThermo's `Jet-A(g)` fit, but keeps
  the C₁₄H₃₀ formation enthalpy and molecular weight, so the burner energy
  balance matches the stoichiometry from [`engine.gasfuel`](@ref) and [`engine.gaschem`](@ref).
- **Agreement with the earlier tables.** Over 200–2500 K, ``c_p`` agrees within
  0.12% for N₂, O₂, CO₂ and Ar, 0.6% for H₂, and 1.8% for H₂O (the old H₂O table
  ran low above 1500 K). Fuel heating values agree within 0.1%. The tests in
  `test/unit_test_thermo_equivalence.jl` check these bounds against a frozen
  copy of the old code in `test/legacy_thermo/`.

### Precomputed mixtures

[`engine.gassum`](@ref) loops over the constituents of a mixture and evaluates
each one. Because the properties are linear in the NASA-9 coefficients, the
mass-fraction-weighted sum of the coefficients describes the mixture exactly.
[`engine.GasMixture`](@ref) holds that sum, and `gassum`, `gassumd` and every
process function built on them accept it in place of the mass-fraction vector
(with the same `n`). Building one costs about as much as one `gassum` call, and
each later property call is one polynomial evaluation. `tfoper!` and `tfsize!`
build one for air once per call, and one each for the combustion products and
the products mixed with cooling air whenever those compositions change.

A `GasMixture` has no entropy of mixing: ``s`` is the mass-weighted sum of the
constituent values, as in `gassum`. An `IdealGasThermo.FrozenGas` built from
mole fractions would add ``-R\sum X_i \ln X_i`` (about 169 J/kg/K for air),
which differs between air and combustion products.

### Adding a species

1. Find the species name in IdealGasThermo's database
   (`IdealGasThermo.spdict.name`). If it is not there, its NASA-9 coefficients
   have to be added to IdealGasThermo first.
2. In `src/engine/gasfun.jl`, add a constant
   `const _IGT_X = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("X"))`,
   pick an unused `igas` index, add a branch for it in [`engine.gasfund`](@ref),
   and add a one-line `gas_X(t)` wrapper like the existing ones.
3. If the species is a fuel, add its C, H, O, N atom counts to
   [`engine.gaschem`](@ref).
4. Add the species to the tests in `test/unit_test_thermo_equivalence.jl`
   (for a fuel, to the heating-value checks).

## Ideal gas with variable $c_p$
```@docs

engine.gas_tset(alpha, n, hspec, tguess)

engine.gas_tsetd(alpha, n, hspec, tguess)

engine.gasfun

engine.gasfund

engine.gaschem(igas)

engine.gassum(alpha, n, t)

engine.gassumd(alpha, n, t)

engine.GasMixture

engine.gas_prat(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)

engine.gas_pratd(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)

engine.gas_pressure_ratio_partials

engine.gas_pressure_ratio_composition_direction

engine.gas_pressure_ratio_composition_directions

engine.gas_delh(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

engine.gas_delhd(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

engine.gas_enthalpy_change_partials

engine.gas_enthalpy_change_composition_direction

engine.gas_enthalpy_change_composition_directions

engine.gas_burn(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

engine.gas_burnd(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

engine.gas_mach(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

engine.gas_machd(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

engine.gas_mass(alpha, n, po, to, ho, so, cpo, ro, mflux, Mguess)

engine.gasfuel(ifuel, n)

engine.gasPr(gas, T)



```
