"""
Thermodynamic wrapper functions that operate on `FlowStation` objects,
delegating to the low-level `gascalc.jl` kernels.

These wrappers form the bridge between the typed-state API (`FlowStation`,
`EngineState`) and the legacy scalar-argument thermodynamics functions
(`gassum`, `gas_prat`, `gas_delh`, `gas_mach`).

## Convention

All mutating wrappers follow the Julia convention of a trailing `!` and
modify their first argument (the station being updated).  They return the
modified station for chaining convenience.

## Species count

`GasState` fixes the species count at 5 (N₂, O₂, CO₂, H₂O, Ar) via
`SVector{5,T}`.  All wrappers pass `n = 5` to the underlying kernels;
this is consistent with the OQ-5 decision recorded in `gas_state.jl`.
"""

const _NAIR = 5   # fixed species count — see OQ-5 in gas_state.jl

# ---------------------------------------------------------------------------
# set_total_from_Tt!
# ---------------------------------------------------------------------------

"""
    set_total_from_Tt!(station::FlowStation) → station

Update the total-state thermodynamic fields of `station` from its current
total temperature `station.Tt` and species composition `station.alpha`.

Specifically, calls `gassum(alpha, 5, Tt)` and stores the results:

| Updated field | From `gassum` output | Meaning                          |
|:--------------|:---------------------|:---------------------------------|
| `station.ht`  | `h`                  | Total complete enthalpy [J/kg]   |
| `station.st`  | `s`                  | Entropy-complement s[Tt] [J/(kg·K)] |
| `station.cpt` | `cp`                 | Specific heat at Tt [J/(kg·K)]   |
| `station.Rt`  | `r`                  | Gas constant [J/(kg·K)]          |

`station.Tt`, `station.pt`, and `station.alpha` are read but **not** modified.
Static-state fields (`Ts`, `ps`, …, `u`) are **not** touched.

## Example

```julia
fs = FlowStation()
fs.alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
fs.Tt = 288.15
set_total_from_Tt!(fs)
# fs.ht, fs.st, fs.cpt, fs.Rt now reflect the thermally-perfect gas at 288.15 K
```

## Note on `st` semantics

`station.st` stores the entropy-complement function `s[Tt]` evaluated at
stagnation temperature.  This is distinct from `station.ss` (= `s[Ts]` at
static temperature), which is populated by `set_static_from_M!`.  Both are
needed for pressure-ratio and Mach-number calculations.
"""
function set_total_from_Tt!(station::FlowStation{T}) where {T<:AbstractFloat}
    s, _, h, _, cp, r = gassum(station.alpha, _NAIR, station.Tt)
    station.ht  = h
    station.st  = s
    station.cpt = cp
    station.Rt  = r
    return station
end

# ---------------------------------------------------------------------------
# set_static_from_M!
# ---------------------------------------------------------------------------

"""
    set_static_from_M!(station::FlowStation, M; epol=1.0) → station

Update the static-state thermodynamic fields of `station` for the given Mach
number `M`, starting from the current total (stagnation) state.

Delegates to `gas_mach(alpha, 5, pt, Tt, ht, st, cpt, Rt, 0.0, M, epol)`,
treating the current station as the stagnation inlet (`mo = 0`).

| Updated field  | Meaning                                      |
|:---------------|:---------------------------------------------|
| `station.Ts`   | Static temperature [K]                       |
| `station.ps`   | Static pressure [Pa]                         |
| `station.hs`   | Static enthalpy [J/kg]                       |
| `station.ss`   | Entropy-complement s[Ts] [J/(kg·K)]          |
| `station.cps`  | Specific heat at static temperature [J/(kg·K)]|
| `station.Rs`   | Gas constant at static conditions [J/(kg·K)] |
| `station.u`    | Flow velocity [m/s]                          |

`station.Tt`, `station.pt`, `station.ht`, `station.st`, `station.cpt`,
`station.Rt`, and `station.alpha` are read but **not** modified.

## Energy conservation

After the call, stagnation enthalpy is conserved:

    station.hs + 0.5 * station.u^2 ≈ station.ht

## Isentropic special case

With `epol = 1.0` (default), the process is isentropic: `station.ss ≈ station.st`.

## Example

```julia
fs = FlowStation()
fs.alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
fs.Tt = 300.0; fs.pt = 1.0e5
set_total_from_Tt!(fs)
set_static_from_M!(fs, 0.8)
# fs.Ts, fs.ps, fs.hs, fs.ss, fs.cps, fs.Rs, fs.u now reflect M=0.8 static state
```
"""
function set_static_from_M!(station::FlowStation{T}, M::Real;
                             epol::Real = one(T)) where {T<:AbstractFloat}
    ps, Ts, hs, ss, cps, Rs = gas_mach(
        station.alpha, _NAIR,
        station.pt, station.Tt, station.ht, station.st, station.cpt, station.Rt,
        zero(T), T(M), T(epol),
    )
    station.Ts  = Ts
    station.ps  = ps
    station.hs  = hs
    station.ss  = ss
    station.cps = cps
    station.Rs  = Rs
    # velocity from Mach: u = M * a = M * sqrt(γ_eff * Rs * Ts)
    # where a² = cps * Rs / (cps - Rs) * Ts  (thermally-perfect gas speed of sound)
    station.u   = T(M) * sqrt(cps * Rs / (cps - Rs) * Ts)
    return station
end
