"""
Nozzle component type for a TASOPT turbofan.

Models the fan nozzle (stations 21→7→8) and the core (turbine exhaust) nozzle
(stations 49c→5→6).  Each nozzle applies a total-pressure loss fraction, then
expands the gas isentropically from the post-loss total state to the ambient
static pressure, detecting choking automatically.

The Nozzle component provides three callable methods:

  - `nozzle_exit`               — forward exit-state computation (no explicit derivatives)
  - `nozzle_massflow_residual`  — continuity residual + partial derivatives for Newton assembly
  - `nozzle_gross_thrust`       — gross momentum + pressure thrust
"""

# ---------------------------------------------------------------------------
# Nozzle struct
# ---------------------------------------------------------------------------

"""
    Nozzle{T<:AbstractFloat}

Parametric nozzle component holding a total-pressure recovery fraction and a
throat area.

Instantiate one `Nozzle` per stream:

```julia
fan_nozzle  = Nozzle(pifn, A7)   # fan nozzle  — uses air   species (alpha)
core_nozzle = Nozzle(pitn, A5)   # core nozzle — uses exhaust species (lambdap)
```

## Fields

| Field | Unit | Description                                                        |
|:------|:-----|:-------------------------------------------------------------------|
| `pn`  | —    | Total-pressure recovery  `pt_nozzle = pn × pt_in`,  0 < pn ≤ 1   |
| `A`   | m²   | Nozzle throat (or exit) area                                       |

## Physics summary

**Step 1 — pressure loss (adiabatic):**

    pt_nozzle = pn × pt_in,    Tt_nozzle = Tt_in,   ht_nozzle = ht_in

**Step 2 — isentropic expansion to ambient:**

If the ideal exit Mach number `M_e < 1` (unchoked):
  exit pressure `p_e = p0`;  `gas_pratd` with `pfn = p0 / pt_nozzle`

If `M_e ≥ 1` (choked):
  exit is sonic, `M_e = 1`;  `gas_machd` with `M = 1`

**Step 3 — exit velocity and density:**

    u_e   = √(2·(ht − h_e))
    ρ_e   = p_e / (R_e · T_e)

**Mass-flow consistency residual** (Newton row for nozzle area matching):

    r = ṁ − ρ_e · u_e · A

**Gross thrust** (momentum + pressure):

    F_gross = ṁ · u_e + (p_e − p0) · A
"""
struct Nozzle{T<:AbstractFloat}
    pn::T   # total-pressure recovery  pt_nozzle = pn × pt_in  [—]
    A ::T   # nozzle throat/exit area                           [m²]
end

# ---------------------------------------------------------------------------
# nozzle_exit — forward exit-state computation
# ---------------------------------------------------------------------------

"""
    nozzle_exit(nozzle, alpha, nair, pt_in, Tt, ht, st, cpt, Rt, p0)
      -> (p_e, T_e, h_e, s_e, cp_e, R_e, u_e, rho_e, M_e)

Compute the nozzle exit static state and exit velocity from the inlet total
conditions, applying the nozzle pressure loss and isentropic expansion.

This is a **forward-mode** function that returns only state values, with no
explicit partial derivatives.  For Newton Jacobian assembly use
`nozzle_massflow_residual` instead.

## Arguments

| Argument | Unit  | Description                                          |
|:---------|:------|:-----------------------------------------------------|
| `nozzle` | —     | `Nozzle` component holding `pn` and `A`              |
| `alpha`  | —     | Species mass-fraction vector (air or exhaust gas)    |
| `nair`   | —     | Number of species (must be 5)                        |
| `pt_in`  | Pa    | Total pressure at the nozzle inlet (pre-loss)        |
| `Tt`     | K     | Total temperature at nozzle inlet                    |
| `ht`     | J/kg  | Total enthalpy at nozzle inlet                       |
| `st`     | J/(kg·K) | Entropy-complement `s[Tt]` at nozzle inlet        |
| `cpt`    | J/(kg·K) | Specific heat at nozzle inlet total state         |
| `Rt`     | J/(kg·K) | Gas constant                                      |
| `p0`     | Pa    | Ambient static pressure                              |

## Returns

`(p_e, T_e, h_e, s_e, cp_e, R_e, u_e, rho_e, M_e)`:

| Return   | Unit     | Description                                        |
|:---------|:---------|:---------------------------------------------------|
| `p_e`    | Pa       | Exit static pressure                               |
| `T_e`    | K        | Exit static temperature                            |
| `h_e`    | J/kg     | Exit static enthalpy                               |
| `s_e`    | J/(kg·K) | Exit entropy-complement                            |
| `cp_e`   | J/(kg·K) | Exit specific heat                                 |
| `R_e`    | J/(kg·K) | Exit gas constant                                  |
| `u_e`    | m/s      | Exit velocity `= √(2·(ht − h_e))`                 |
| `rho_e`  | kg/m³    | Exit density `= p_e / (R_e · T_e)`                |
| `M_e`    | —        | Exit Mach number (capped at 1 if choked)           |
"""
function nozzle_exit(
    nozzle::Nozzle,
    alpha, nair::Int,
    pt_in ::T,
    Tt    ::T,
    ht    ::T,
    st    ::T,
    cpt   ::T,
    Rt    ::T,
    p0    ::T,
) where {T<:AbstractFloat}
    # gas_pratd / gas_machd use Float64 arithmetic internally (zeros(n), Newton iteration).
    # Convert to Float64 for the gas function calls; results are converted back to T.
    pt_d  = Float64(nozzle.pn) * Float64(pt_in)
    Tt_d  = Float64(Tt)
    ht_d  = Float64(ht)
    st_d  = Float64(st)
    cpt_d = Float64(cpt)
    Rt_d  = Float64(Rt)
    p0_d  = Float64(p0)

    # Step 2 — probe unchoked expansion to p0
    pfn = p0_d / pt_d
    p_e64, T_e64, h_e64, s_e64, cp_e64, R_e64,
    _p_po,
    _p_so, _T_so, _h_so, _s_so,
    _p_pi, _T_pi, _h_pi, _s_pi,
    _p_ep, _T_ep, _h_ep, _s_ep,
    _p_al, _T_al, _h_al, _s_al, _cp_al, _R_al = gas_pratd(alpha, nair,
        pt_d, Tt_d, ht_d, st_d, cpt_d, Rt_d, pfn, 1.0)

    u_e64 = sqrt(2 * max(ht_d - h_e64, 0.0))
    M_e64 = u_e64 / sqrt(T_e64 * cp_e64 * R_e64 / (cp_e64 - R_e64))

    if M_e64 >= 1.0
        # Choked: re-solve at M = 1
        p_e64, T_e64, h_e64, s_e64, cp_e64, R_e64,
        _p_so, _p_po, _p_ep,
        _p_Tt, _T_Tt, _h_Tt, _s_Tt,
        _p_ht, _T_ht, _h_ht, _s_ht,
        _p_M, _T_M, _h_M, _s_M,
        _p_al, _T_al, _h_al, _s_al, _cp_al, _R_al = gas_machd(alpha, nair,
            pt_d, Tt_d, ht_d, st_d, cpt_d, Rt_d, 0.0, 1.0, 1.0)

        u_e64 = sqrt(2 * max(ht_d - h_e64, 0.0))
        M_e64 = 1.0
    end

    rho_e64 = p_e64 / (R_e64 * T_e64)
    # Convert results back to T
    return T(p_e64), T(T_e64), T(h_e64), T(s_e64), T(cp_e64), T(R_e64),
           T(u_e64), T(rho_e64), T(M_e64)
end

# ---------------------------------------------------------------------------
# nozzle_massflow_residual — continuity residual + derivatives for Newton
# ---------------------------------------------------------------------------

"""
    nozzle_massflow_residual(nozzle, alpha, nair, pt_in, Tt, ht, st, cpt, Rt, p0, mdot)
      -> (r, r_pt, r_ht, r_st, r_Tt, p_al, T_al, h_al, s_al, cp_al, R_al)

Compute the nozzle continuity residual

    r = ṁ − ρ_e · u_e · A

and its partial derivatives with respect to the nozzle-inlet total-state
variables `(pt_in, ht, st, Tt)` for Newton Jacobian assembly.

## Arguments

| Argument | Unit     | Description                                         |
|:---------|:---------|:----------------------------------------------------|
| `nozzle` | —        | `Nozzle` component                                  |
| `alpha`  | —        | Species mass-fraction vector                        |
| `nair`   | —        | Number of species (must be 5)                       |
| `pt_in`  | Pa       | Total pressure at the nozzle inlet (pre-loss)       |
| `Tt`     | K        | Total temperature at nozzle inlet                   |
| `ht`     | J/kg     | Total enthalpy at nozzle inlet                      |
| `st`     | J/(kg·K) | Entropy-complement at nozzle inlet                  |
| `cpt`    | J/(kg·K) | Specific heat at nozzle inlet total state           |
| `Rt`     | J/(kg·K) | Gas constant                                        |
| `p0`     | Pa       | Ambient static pressure                             |
| `mdot`   | kg/s     | Inlet mass flow rate                                |

## Returns

`(r, r_pt, r_ht, r_st, r_Tt, p_al, T_al, h_al, s_al, cp_al, R_al)`:

| Return | Description                                                         |
|:-------|:--------------------------------------------------------------------|
| `r`    | Continuity residual `= ṁ − ρ_e · u_e · A`                         |
| `r_pt` | `∂r/∂pt_in = −A · ∂(ρ_e u_e)/∂pt_in`                             |
| `r_ht` | `∂r/∂ht    = −A · ∂(ρ_e u_e)/∂ht`                                |
| `r_st` | `∂r/∂st    = −A · ∂(ρ_e u_e)/∂st`                                |
| `r_Tt` | `∂r/∂Tt    = −A · ∂(ρ_e u_e)/∂Tt` (non-zero only when choked)    |
| `p_al` | Species sensitivities `∂p_e/∂αᵢ` for fuel-fraction chain rule     |
| `T_al` | Species sensitivities `∂T_e/∂αᵢ`                                  |
| `h_al` | Species sensitivities `∂h_e/∂αᵢ`                                  |
| `s_al` | Species sensitivities `∂s_e/∂αᵢ`                                  |
| `cp_al`| Species sensitivities `∂cp_e/∂αᵢ`                                 |
| `R_al` | Species sensitivities `∂R_e/∂αᵢ`                                  |

## Newton Jacobian chain rule

The Newton driver assembles the Jacobian row for the nozzle continuity equation
by chaining the returned partial derivatives with the upstream sensitivities
`(pt_in_x, ht_x, st_x, Tt_x)` for each Newton variable `x`:

```
a[row, x] = r_pt · pt_in_x + r_ht · ht_x + r_st · st_x + r_Tt · Tt_x + ...
```

plus species contributions via `p_al`, `T_al`, `h_al`, `s_al` for the
fuel-fraction chain rule (fuel-to-air ratio `ff`, offtake fraction `fo`).

## Choked vs. unchoked derivative structure

**Unchoked** (`M_e < 1`):  exit state depends on `(pt_in, st)` only;
`r_Tt = 0` exactly.  `ht` enters only through the energy equation `u = √(2·(ht−h_e))`,
so `r_ht = −ρ_e · A / u_e`.

**Choked** (`M_e = 1`):  the sonic condition couples the exit state to all four
total-state variables `(pt_in, ht, Tt, st)`, so all four partial derivatives
are generally nonzero.
"""
function nozzle_massflow_residual(
    nozzle::Nozzle,
    alpha, nair::Int,
    pt_in ::T,
    Tt    ::T,
    ht    ::T,
    st    ::T,
    cpt   ::T,
    Rt    ::T,
    p0    ::T,
    mdot  ::T,
) where {T<:AbstractFloat}
    A  = nozzle.A
    pn = nozzle.pn

    # gas_pratd / gas_machd require Float64 arithmetic internally.
    # Convert all scalar state inputs to Float64; derivatives are computed in
    # Float64 and converted back to T for the returned residual partials.
    pt64 = Float64(pn) * Float64(pt_in)   # post-loss total pressure
    Tt64 = Float64(Tt)
    ht64 = Float64(ht)
    st64 = Float64(st)
    cpt64 = Float64(cpt)
    Rt64  = Float64(Rt)
    p064  = Float64(p0)

    # ---------- probe unchoked (p_e = p0) expansion ----------
    pfn    = p064 / pt64
    pfn_pt = -pfn / pt64       # ∂pfn/∂pt_nozzle  (Float64)

    p_e, T_e, h_e, s_e, cp_e, R_e,
    p_e_pt_d,                        # ∂p_e/∂pt at fixed pfn, st
    p_e_st, T_e_st, h_e_st, s_e_st,  # ∂/∂st
    p_e_pfn, T_e_pfn, h_e_pfn, s_e_pfn,  # ∂/∂pfn
    _p_ep, _T_ep, _h_ep, _s_ep,      # ∂/∂epol (not needed, epol=1)
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(alpha, nair,
        pt64, Tt64, ht64, st64, cpt64, Rt64, pfn, 1.0)

    u_e0 = sqrt(2 * max(ht64 - h_e, 0.0))
    M_e  = u_e0 / sqrt(T_e * cp_e * R_e / (cp_e - R_e))

    # Declare before if-else to ensure they are in scope after the branch (Float64 working vars)
    rhou_pt_nozzle = 0.0   # ∂(ρ_e u_e)/∂pt_nozzle; assigned in each branch
    rhou_ht        = 0.0   # ∂(ρ_e u_e)/∂ht
    rhou_st        = 0.0   # ∂(ρ_e u_e)/∂st
    rhou_Tt        = 0.0   # ∂(ρ_e u_e)/∂Tt
    rho_e          = 0.0   # exit density
    u_e            = 0.0   # exit velocity

    if M_e < 1.0
        # ---- Unchoked ----
        u_e  = u_e0

        # Total ∂p/∂pt_nozzle (add pfn chain to direct derivative)
        p_e_pt = p_e_pfn * pfn_pt + p_e_pt_d
        T_e_pt = T_e_pfn * pfn_pt
        h_e_pt = h_e_pfn * pfn_pt
        # h_e_pt_d ≈ 0 for gas_pratd at epol=1 (p=po*pfn, independent of po alone)

        rho_e = p_e / (R_e * T_e)

        # Guard against u_e → 0 (static or near-static conditions)
        u_safe = max(u_e, 0.02 * sqrt(Rt64 * Tt64))

        # Unchoked: h_e does NOT depend on ht (gas_pratd only uses so, pratio, po)
        # so u_e_ht = 1/u_e (direct ht contribution to energy equation only)
        # rho_e does NOT depend on ht either.
        rhou_pt_nozzle = rho_e * (-h_e_pt / u_safe) +
                         u_e   * (p_e_pt / (R_e * T_e) - rho_e / T_e * T_e_pt)
        rhou_ht        = rho_e / u_safe
        rhou_st        = rho_e * (-h_e_st / u_safe) +
                         u_e   * (p_e_st / (R_e * T_e) - rho_e / T_e * T_e_st)
        rhou_Tt        = 0.0

    else
        # ---- Choked (M = 1) ----
        p_e, T_e, h_e, s_e, cp_e, R_e,
        p_e_st,              # ∂p_e/∂st
        p_e_pt,              # ∂p_e/∂pt_nozzle
        _p_ep,
        p_e_Tt, T_e_Tt, h_e_Tt, s_e_Tt,  # ∂/∂Tt
        p_e_ht, T_e_ht, h_e_ht, s_e_ht,  # ∂/∂ht
        _p_M, _T_M, _h_M, _s_M,
        p_al, T_al, h_al, s_al,
        cp_al, R_al = gas_machd(alpha, nair,
            pt64, Tt64, ht64, st64, cpt64, Rt64, 0.0, 1.0, 1.0)

        u_e   = sqrt(2 * max(ht64 - h_e, 0.0))
        rho_e = p_e / (R_e * T_e)

        u_safe = max(u_e, 0.02 * sqrt(Rt64 * Tt64))

        # Choked: T_e and h_e depend on (ht, Tt); p_e also depends on (pt, st)
        # In choked case T_e does NOT depend on (pt_nozzle, st) directly.
        rhou_pt_nozzle = u_e * p_e_pt / (R_e * T_e)
        rhou_st        = u_e * p_e_st / (R_e * T_e)
        rhou_ht        = u_e * (p_e_ht / (R_e * T_e) - rho_e / T_e * T_e_ht) +
                         rho_e * (1.0 - h_e_ht) / u_safe
        rhou_Tt        = u_e * (p_e_Tt / (R_e * T_e) - rho_e / T_e * T_e_Tt) +
                         rho_e * (-h_e_Tt) / u_safe
    end

    # Convert ∂/∂pt_nozzle → ∂/∂pt_in  (pt_nozzle = pn × pt_in)
    # All derivative computations above are in Float64; convert to T
    pn_f64 = Float64(pn)
    rhou_pt = T(rhou_pt_nozzle * pn_f64)
    rhou_ht = T(rhou_ht)
    rhou_st = T(rhou_st)
    rhou_Tt = T(rhou_Tt)
    rho_e_T = T(rho_e)
    u_e_T   = T(u_e)

    # Residual and its partial derivatives
    rhoUA = rho_e_T * u_e_T * A
    r     = T(Float64(mdot) - Float64(rho_e) * Float64(u_e) * Float64(A))

    r_pt = -rhou_pt * A
    r_ht = -rhou_ht * A
    r_st = -rhou_st * A
    r_Tt = -rhou_Tt * A

    return r, r_pt, r_ht, r_st, r_Tt, p_al, T_al, h_al, s_al, cp_al, R_al
end

# ---------------------------------------------------------------------------
# nozzle_gross_thrust — momentum + pressure thrust
# ---------------------------------------------------------------------------

"""
    nozzle_gross_thrust(nozzle, mdot, u_e, p_e, p0) -> F_gross

Compute the gross thrust contribution of one nozzle stream:

    F_gross = ṁ · u_e + (p_e − p0) · A

where `A = nozzle.A` is the nozzle throat area.

The net engine thrust is assembled by the caller from the fan and core
gross thrust contributions, minus the ram drag.  For an ideally expanded
nozzle (`p_e = p0`) the pressure term vanishes and `F_gross = ṁ · u_e`.

## Arguments

| Argument | Unit  | Description                          |
|:---------|:------|:-------------------------------------|
| `nozzle` | —     | `Nozzle` component (uses `A` field)  |
| `mdot`   | kg/s  | Mass flow rate through the nozzle    |
| `u_e`    | m/s   | Exit velocity                        |
| `p_e`    | Pa    | Exit static pressure                 |
| `p0`     | Pa    | Ambient static pressure              |

## Returns

`F_gross` [N]: gross thrust from this nozzle stream.
"""
function nozzle_gross_thrust(
    nozzle::Nozzle,
    mdot  ::T,
    u_e   ::T,
    p_e   ::T,
    p0    ::T,
) where {T<:AbstractFloat}
    return mdot * u_e + (p_e - p0) * nozzle.A
end
