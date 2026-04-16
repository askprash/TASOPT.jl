"""
Shaft component types for a TASOPT turbofan.

Models the HP and LP mechanical power transmission shafts, encapsulating the
spool torque (power) balance between turbines and compressors and the gear-ratio
speed constraint between the fan and LPC.
"""

# ---------------------------------------------------------------------------
# Shaft struct
# ---------------------------------------------------------------------------

"""
    Shaft{T<:AbstractFloat}

Parametric shaft component holding the mechanical power loss fraction and the
gear ratio.

Instantiate one `Shaft` per spool:

```julia
shaft_hp = Shaft(epsh, T(1.0))          # HP shaft: HPT drives HPC
shaft_lp = Shaft(epsl, Gearf)           # LP shaft: LPT drives fan + LPC
```

## Fields

| Field   | Unit | Description                                                     |
|:--------|:-----|:----------------------------------------------------------------|
| `eps`   | —    | Mechanical power loss fraction  (0 ≤ eps < 1)                  |
| `Gearf` | —    | Gear ratio Nl/Nf  (= 1.0 for HP shaft, physical ratio for LP)  |

## Physics

The HP shaft power balance (HPT drives HPC):

    W_HPT · (1 − eps) = W_HPC

The LP shaft power balance (LPT drives fan + LPC + offtake):

    W_LPT · (1 − eps) = W_fan + W_LPC + W_offtake

The LP shaft speed constraint (gear-ratio link between fan and LPC):

    Gearf · trf · Nf − trl · Nl = 0
"""
struct Shaft{T<:AbstractFloat}
    eps  ::T    # mechanical power loss fraction  [—]
    Gearf::T    # gear ratio Nl/Nf               [—]
end

# ---------------------------------------------------------------------------
# hp_shaft_work — HPT specific work from HPC demand
# ---------------------------------------------------------------------------

"""
    hp_shaft_work(shaft, fo, ff, ht3, ht25c)
      -> (dhht, dhfac, dhfac_fo, dhfac_ff)

Compute the HPT specific work required to drive the HPC on the HP shaft,
returning the value and the key scaling factors needed to assemble the Newton
Jacobian by chain rule in the calling solver.

## Physics

The HP shaft power balance imposes:

    dhfac = −(1 − fo) / (1 − fo + ff) / (1 − shaft.eps)
    dhht  = (ht3 − ht25c) · dhfac

where `fo` is the normalised offtake mass-flow fraction, `ff` the fuel/air
ratio, `ht3` the HPC exit enthalpy, and `ht25c` the HPC inlet (intercooler
exit) enthalpy.

## Arguments

| Argument | Unit  | Description                                     |
|:---------|:------|:------------------------------------------------|
| `shaft`  | —     | `Shaft` component (HP shaft with `Gearf = 1`)   |
| `fo`     | —     | Normalised bleed/offtake mass-flow fraction     |
| `ff`     | —     | Fuel/air ratio                                  |
| `ht3`    | J/kg  | HPC exit total enthalpy                         |
| `ht25c`  | J/kg  | HPC inlet (post-intercooler) total enthalpy     |

## Returns

`(dhht, dhfac, dhfac_fo, dhfac_ff)`:

| Return     | Description                                    |
|:-----------|:-----------------------------------------------|
| `dhht`     | HPT specific work  (< 0, enthalpy drop)        |
| `dhfac`    | Scaling factor  = dhht / (ht3 − ht25c)         |
| `dhfac_fo` | ∂dhfac/∂fo                                     |
| `dhfac_ff` | ∂dhfac/∂ff                                     |

The Newton Jacobian entry `dhht_x` for any Newton variable `x` is assembled
by the caller as:

    dhht_x = (ht3 − ht25c) · dhfac_x  +  dhfac · (ht3_x − ht25c_x)

where `dhfac_x = dhfac_fo · fo_x + dhfac_ff · ff_x`.

## Design-point identity

At the design operating point, `dhht` exactly equals the HPC work divided
by the mechanical efficiency `(1 − eps)`.
"""
function hp_shaft_work(
    shaft ::Shaft{T},
    fo    ::T,
    ff    ::T,
    ht3   ::T,
    ht25c ::T,
) where {T<:AbstractFloat}
    fac      = one(T) - fo + ff
    dhfac    = -(one(T) - fo) / fac / (one(T) - shaft.eps)
    dhfac_fo =  dhfac / fac + one(T) / fac / (one(T) - shaft.eps)
    dhfac_ff = -dhfac / fac

    dhht = (ht3 - ht25c) * dhfac

    return dhht, dhfac, dhfac_fo, dhfac_ff
end

# ---------------------------------------------------------------------------
# lp_shaft_work — LPT specific work from fan + LPC demand
# ---------------------------------------------------------------------------

"""
    lp_shaft_work(shaft, fo, ff, BPR, ht25, ht19c, ht21, ht2, Pom)
      -> (dhlt, dlfac, dlfac_fo, dlfac_ff)

Compute the LPT specific work required to drive the fan, LPC, and power
offtake on the LP shaft, returning the value and the key scaling factors
needed to assemble the Newton Jacobian by chain rule in the calling solver.

## Physics

The LP shaft power balance imposes:

    dlfac  = −1 / (1 − fo + ff) / (1 − shaft.eps)
    demand = ht25 − ht19c + BPR · (ht21 − ht2) + Pom
    dhlt   = demand · dlfac

where `ht25` and `ht19c` are the LPC exit and inlet enthalpies, `ht21` and
`ht2` are the fan exit and inlet enthalpies, and `Pom` is the normalised
power offtake per unit core mass flow.

## Arguments

| Argument | Unit  | Description                                          |
|:---------|:------|:-----------------------------------------------------|
| `shaft`  | —     | `Shaft` component (LP shaft with physical `Gearf`)   |
| `fo`     | —     | Normalised bleed/offtake mass-flow fraction          |
| `ff`     | —     | Fuel/air ratio                                       |
| `BPR`    | —     | Bypass ratio  mf / mcore                             |
| `ht25`   | J/kg  | LPC exit total enthalpy                              |
| `ht19c`  | J/kg  | LPC inlet (post-precooler) total enthalpy            |
| `ht21`   | J/kg  | Fan exit total enthalpy                              |
| `ht2`    | J/kg  | Fan inlet (engine face) total enthalpy               |
| `Pom`    | J/kg  | Normalised power offtake per unit core mass flow     |

## Returns

`(dhlt, dlfac, dlfac_fo, dlfac_ff)`:

| Return     | Description                                         |
|:-----------|:----------------------------------------------------|
| `dhlt`     | LPT specific work  (< 0, enthalpy drop)             |
| `dlfac`    | Scaling factor  = dhlt / demand                     |
| `dlfac_fo` | ∂dlfac/∂fo                                          |
| `dlfac_ff` | ∂dlfac/∂ff                                          |

The Newton Jacobian entries are assembled by the caller using these factors
together with the demand-component sensitivities (∂ht25/∂x, BPR_x, etc.).
"""
function lp_shaft_work(
    shaft ::Shaft{T},
    fo    ::T,
    ff    ::T,
    BPR   ::T,
    ht25  ::T,
    ht19c ::T,
    ht21  ::T,
    ht2   ::T,
    Pom   ::T,
) where {T<:AbstractFloat}
    fac      = one(T) - fo + ff
    dlfac    = -one(T) / fac / (one(T) - shaft.eps)
    dlfac_fo =  dlfac / fac
    dlfac_ff = -dlfac / fac

    demand = ht25 - ht19c + BPR * (ht21 - ht2) + Pom
    dhlt   = demand * dlfac

    return dhlt, dlfac, dlfac_fo, dlfac_ff
end

# ---------------------------------------------------------------------------
# shaft_speed_residual — gear-ratio speed constraint for the LP shaft
# ---------------------------------------------------------------------------

"""
    shaft_speed_residual(shaft, Nf, Nl, trf, trl)
      -> (r, r_Nf, r_Nl)

Speed-constraint residual for the LP shaft gear-ratio link between the fan
(Nf) and the LPC (Nl) corrected speeds.

## Physics

The gear ratio fixes the ratio of LPC-spool speed to fan speed:

    r = shaft.Gearf · trf · Nf − trl · Nl = 0

where `trf = √(Tt2 / Tref)` and `trl = √(Tt19c / Tref)` are the temperature
correction factors for the fan and LPC inlet stations respectively.

## Arguments

| Argument | Unit | Description                                       |
|:---------|:-----|:--------------------------------------------------|
| `shaft`  | —    | `Shaft` component (LP shaft with physical `Gearf`)|
| `Nf`     | —    | Fan corrected spool speed                         |
| `Nl`     | —    | LPC corrected spool speed                         |
| `trf`    | —    | Fan inlet temperature correction √(Tt2/Tref)      |
| `trl`    | —    | LPC inlet temperature correction √(Tt19c/Tref)    |

## Returns

`(r, r_Nf, r_Nl)` — the residual and its partial derivatives w.r.t. `Nf`
and `Nl`, for use by the Newton Jacobian assembly (row 1 of the 9×9 system).

## Design-point identity

At the design operating point, `Gearf · trf · Nf = trl · Nl` exactly, so
`r = 0`.
"""
# ---------------------------------------------------------------------------
# hp_shaft_workd — HPT work with product-rule Jacobian pre-applied
# ---------------------------------------------------------------------------

"""
    hp_shaft_workd(shaft, fo, ff, ht3, ht25c)
        -> (dhht, dhht_fo, dhht_ff, dhht_ht3, dhht_ht25c)

Compute the HPT specific work and its partial derivatives with respect to
the four natural inputs `(fo, ff, ht3, ht25c)`.

Wraps `hp_shaft_work` and applies the product rule for
`dhht = (ht3 − ht25c) · dhfac(fo, ff)` so that the caller assembles Newton
Jacobian entries by simple chain rule:

    dhht_x = dhht_fo · fo_x  +  dhht_ff · ff_x
           + dhht_ht3 · ht3_x  +  dhht_ht25c · ht25c_x

## Returns

| Return       | Description                        |
|:-------------|:-----------------------------------|
| `dhht`       | HPT specific work  (< 0)          |
| `dhht_fo`    | ∂dhht/∂fo                          |
| `dhht_ff`    | ∂dhht/∂ff                          |
| `dhht_ht3`   | ∂dhht/∂ht3   = dhfac               |
| `dhht_ht25c` | ∂dhht/∂ht25c = −dhfac              |
"""
function hp_shaft_workd(
    shaft ::Shaft{T},
    fo    ::T,
    ff    ::T,
    ht3   ::T,
    ht25c ::T,
) where {T<:AbstractFloat}
    dhht, dhfac, dhfac_fo, dhfac_ff = hp_shaft_work(shaft, fo, ff, ht3, ht25c)
    Δh = ht3 - ht25c
    dhht_fo    = Δh * dhfac_fo
    dhht_ff    = Δh * dhfac_ff
    dhht_ht3   = dhfac
    dhht_ht25c = -dhfac
    return dhht, dhht_fo, dhht_ff, dhht_ht3, dhht_ht25c
end

# ---------------------------------------------------------------------------
# lp_shaft_workd — LPT work with product-rule Jacobian pre-applied
# ---------------------------------------------------------------------------

"""
    lp_shaft_workd(shaft, fo, ff, BPR, ht25, ht19c, ht21, ht2, Pom)
        -> (dhlt, dhlt_fo, dhlt_ff, dhlt_BPR, dhlt_ht25, dhlt_ht19c,
            dhlt_ht21, dhlt_ht2, dhlt_Pom)

Compute the LPT specific work and its partial derivatives with respect to
the nine natural inputs.

Wraps `lp_shaft_work` and applies the product rule for
`dhlt = demand · dlfac(fo, ff)` where
`demand = (ht25 − ht19c) + BPR · (ht21 − ht2) + Pom`,
so the caller assembles Newton entries by chain rule:

    dhlt_x = dhlt_fo · fo_x  +  dhlt_ff · ff_x  +  dhlt_BPR · BPR_x
           + dhlt_ht25 · ht25_x  +  dhlt_ht19c · ht19c_x
           + dhlt_ht21 · ht21_x  +  dhlt_ht2 · ht2_x  +  dhlt_Pom · Pom_x

## Returns

| Return       | Description                                     |
|:-------------|:------------------------------------------------|
| `dhlt`       | LPT specific work  (< 0)                       |
| `dhlt_fo`    | ∂dhlt/∂fo   = demand · dlfac_fo                 |
| `dhlt_ff`    | ∂dhlt/∂ff   = demand · dlfac_ff                 |
| `dhlt_BPR`   | ∂dhlt/∂BPR  = (ht21 − ht2) · dlfac             |
| `dhlt_ht25`  | ∂dhlt/∂ht25 = dlfac                             |
| `dhlt_ht19c` | ∂dhlt/∂ht19c = −dlfac                           |
| `dhlt_ht21`  | ∂dhlt/∂ht21 = BPR · dlfac                       |
| `dhlt_ht2`   | ∂dhlt/∂ht2  = −BPR · dlfac                      |
| `dhlt_Pom`   | ∂dhlt/∂Pom  = dlfac                             |
"""
function lp_shaft_workd(
    shaft ::Shaft{T},
    fo    ::T,
    ff    ::T,
    BPR   ::T,
    ht25  ::T,
    ht19c ::T,
    ht21  ::T,
    ht2   ::T,
    Pom   ::T,
) where {T<:AbstractFloat}
    dhlt, dlfac, dlfac_fo, dlfac_ff = lp_shaft_work(
        shaft, fo, ff, BPR, ht25, ht19c, ht21, ht2, Pom)

    demand     = ht25 - ht19c + BPR * (ht21 - ht2) + Pom
    dhlt_fo    = demand * dlfac_fo
    dhlt_ff    = demand * dlfac_ff
    dhlt_BPR   = (ht21 - ht2) * dlfac
    dhlt_ht25  = dlfac
    dhlt_ht19c = -dlfac
    dhlt_ht21  = BPR * dlfac
    dhlt_ht2   = -BPR * dlfac
    dhlt_Pom   = dlfac

    return dhlt, dhlt_fo, dhlt_ff, dhlt_BPR,
           dhlt_ht25, dhlt_ht19c, dhlt_ht21, dhlt_ht2, dhlt_Pom
end

# ---------------------------------------------------------------------------
# shaft_speed_residual — gear-ratio speed constraint for the LP shaft
# ---------------------------------------------------------------------------

function shaft_speed_residual(
    shaft ::Shaft{T},
    Nf    ::T,
    Nl    ::T,
    trf   ::T,
    trl   ::T,
) where {T<:AbstractFloat}
    r    = shaft.Gearf * trf * Nf - trl * Nl
    r_Nf = shaft.Gearf * trf
    r_Nl = -trl
    return r, r_Nf, r_Nl
end
