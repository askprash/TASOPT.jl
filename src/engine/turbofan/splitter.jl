"""
Splitter component type for a TASOPT turbofan.

Models the fan-outlet junction that divides fan discharge into the bypass
stream (stations 7–8, fan nozzle) and the core stream (station 19c, LPC inlet).
The bypass ratio (BPR) is determined kinematically from the corrected mass
flows at the fan and LPC inlets and the local total-state conditions.
"""

# ---------------------------------------------------------------------------
# Splitter struct
# ---------------------------------------------------------------------------

"""
    Splitter

Stateless splitter component representing the fan-outlet flow junction.

The `Splitter` type owns the bypass-ratio calculation and its partial
derivatives with respect to the Newton solver's independent variables.
It carries no design fields because BPR is a *derived* quantity that
emerges from the corrected mass-flow matching between the fan and LPC —
there is no independent BPR design input to tfoper.

## Physics summary

The bypass ratio is the ratio of bypass corrected mass flow to core corrected
mass flow, both referred to the same reference state `(Tref, pref)`:

```
BPR = (mf / ml) × √(Tt19c / Tt2) × (pt2 / pt19c)
```

where `mf` and `ml` are the dimensionless corrected fan and LPC mass flows,
and `Tt2 / pt2`, `Tt19c / pt19c` are the fan-inlet and LPC-inlet total states.

Partial derivatives with respect to the direct arguments are:

```
∂BPR/∂mf    =  BPR / mf      (explicit mf dependence)
∂BPR/∂ml    = −BPR / ml      (explicit ml dependence)
∂BPR/∂pt2   =  BPR / pt2     (through corrected mass-flow ratio)
∂BPR/∂pt19c = −BPR / pt19c   (through corrected mass-flow ratio)
```

The caller (Newton Jacobian assembly) uses these together with the upstream
sensitivities `pt2_mf`, `pt19c_mf`, `pt2_ml`, `pt19c_ml`, `pt2_Mi`,
`pt19c_Mi` (from the BLI mixing calculation) to assemble the full derivatives:

```
BPR_mf = ∂BPR/∂mf  + BPR_pt2 × pt2_mf  + BPR_pt19c × pt19c_mf
BPR_ml = ∂BPR/∂ml  + BPR_pt2 × pt2_ml  + BPR_pt19c × pt19c_ml
BPR_Mi =             BPR_pt2 × pt2_Mi   + BPR_pt19c × pt19c_Mi
```
"""
struct Splitter end

# ---------------------------------------------------------------------------
# bypass_ratio — kinematic BPR from corrected mass flows
# ---------------------------------------------------------------------------

"""
    bypass_ratio(splitter, mf, ml, pt2, pt19c, Tt2, Tt19c)
      -> (BPR, BPR_mf, BPR_ml, BPR_pt2, BPR_pt19c)

Compute the bypass ratio and its partial derivatives with respect to the
direct arguments, for use in Newton Jacobian assembly.

## Arguments

| Argument  | Unit | Description                                                    |
|:----------|:-----|:---------------------------------------------------------------|
| `splitter`| —    | `Splitter` component (stateless, used for dispatch)            |
| `mf`      | —    | Dimensionless corrected fan mass flow (Newton iterate)         |
| `ml`      | —    | Dimensionless corrected LPC mass flow (Newton iterate)         |
| `pt2`     | Pa   | Fan-inlet total pressure  (function of `mf`, `ml`, `Mi`)      |
| `pt19c`   | Pa   | LPC-inlet total pressure  (function of `mf`, `ml`, `Mi`)      |
| `Tt2`     | K    | Fan-inlet total temperature  (set from freestream; not a Newton iterate) |
| `Tt19c`   | K    | LPC-inlet total temperature  (set from freestream; not a Newton iterate) |

## Returns

`(BPR, BPR_mf, BPR_ml, BPR_pt2, BPR_pt19c)`:

| Return      | Description                                                   |
|:------------|:--------------------------------------------------------------|
| `BPR`       | Bypass ratio  = mf/ml × √(Tt19c/Tt2) × pt2/pt19c            |
| `BPR_mf`    | ∂BPR/∂mf  (direct, without chain through pt2/pt19c)          |
| `BPR_ml`    | ∂BPR/∂ml  (direct, without chain through pt2/pt19c)          |
| `BPR_pt2`   | ∂BPR/∂pt2                                                     |
| `BPR_pt19c` | ∂BPR/∂pt19c                                                   |

The Newton Jacobian entries `BPR_x` for any Newton variable `x` are assembled
by the caller as:

    BPR_x = BPR_mf × mf_x + BPR_ml × ml_x
          + BPR_pt2 × pt2_x + BPR_pt19c × pt19c_x

where `mf_x = 1` (if x = mf, else 0) and `pt2_x`, `pt19c_x` are the upstream
sensitivities from the BLI mixing / inlet calculation.

## Design-point identity

At the design corrected mass flows `(mf_D, ml_D)` and design total states, BPR
equals the design bypass ratio.  The identity `BPR_mf × ml = BPR` holds
algebraically (Euler's homogeneous-function theorem, since BPR is homogeneous
of degree 1 in `mf / ml`).
"""
function bypass_ratio(
    ::Splitter,
    mf    ::T,
    ml    ::T,
    pt2   ::T,
    pt19c ::T,
    Tt2   ::T,
    Tt19c ::T,
) where {T<:AbstractFloat}
    # Common sub-expression: √(Tt19c/Tt2) × pt2/pt19c
    sqr       = sqrt(Tt19c / Tt2) * pt2 / pt19c
    BPR       = mf / ml * sqr
    # BPR_mf is computed as 1/ml × sqr (not BPR/mf) to match the floating-point
    # evaluation order of the original inline code in tfoper.jl and thereby
    # preserve bit-for-bit regression against the reference solution.
    BPR_mf    =  one(T) / ml * sqr
    BPR_ml    = -BPR / ml
    BPR_pt2   =  BPR / pt2
    BPR_pt19c = -BPR / pt19c
    return BPR, BPR_mf, BPR_ml, BPR_pt2, BPR_pt19c
end
