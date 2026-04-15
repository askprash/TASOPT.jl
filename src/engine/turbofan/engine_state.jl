"""
Engine-state container: all named flow stations, ambient flight condition,
and a composed design-state for a TASOPT turbofan.
"""

# ---------------------------------------------------------------------------
# EngineState
# ---------------------------------------------------------------------------

"""
    EngineState{T<:AbstractFloat}

Mutable container for the complete aerothermodynamic state of a TASOPT
turbofan at a single operating point.  It holds:

- **One `FlowStation{T}` per enumerated `EngineStation`**, named by the
  TASOPT station number with an `st` prefix: `st0`, `st2`, `st3`, …,
  `st49c`, etc.  Stations with fractional or suffixed numbers use the
  same string — `st19c`, `st4a`, `st49c` — which are all valid Julia
  identifiers.

- **Ambient / flight-condition scalars** that describe the far-field
  operating point:

  | Field | Unit  | Description                         |
  |:------|:------|:------------------------------------|
  | `M0`  | —     | Freestream Mach number              |
  | `T0`  | K     | Freestream static temperature       |
  | `p0`  | Pa    | Freestream static pressure          |
  | `a0`  | m/s   | Freestream speed of sound           |

- **A composed `DesignState{T}`** field (`design`) carrying the frozen
  design-point scalars that are set once during `tfsize!` and read (but
  never mutated) by `tfoper!`.

## Station fields

| Field    | `EngineStation` member | TASOPT # | Description                        |
|:---------|:-----------------------|:--------:|:-----------------------------------|
| `st0`    | `Freestream`           | 0        | Far-field ambient                  |
| `st18`   | `FanFaceOuter`         | 18       | Fan face, outside casing BL        |
| `st19`   | `FanFaceLPC`           | 19       | Fan face, LPC stream               |
| `st19c`  | `PreCoolerOut`         | 19c      | Pre-cooler outlet / LPC inlet      |
| `st2`    | `FanFaceFan`           | 2        | Fan face, fan stream               |
| `st21`   | `FanExit`              | 21       | Fan exit / pre-cooler inlet        |
| `st25`   | `LPCExit`              | 25       | LPC exit / inter-cooler inlet      |
| `st25c`  | `InterCoolerOut`       | 25c      | Inter-cooler outlet / HPC inlet    |
| `st3`    | `HPCExit`              | 3        | HPC exit / combustor inlet         |
| `st4`    | `CombustorExit`        | 4        | Combustor exit (before cooling)    |
| `st4a`   | `CoolMixInlet`         | 4a       | Cooling-flow outlet / mix inlet    |
| `st41`   | `TurbineInlet`         | 41       | Turbine inlet (after cooling mix)  |
| `st45`   | `HPTExit`              | 45       | HPT exit / LPT inlet               |
| `st49`   | `LPTExit`              | 49       | LPT exit / regen-cooler inlet      |
| `st49c`  | `RegenCoolerOut`       | 49c      | Regenerative cooler outlet         |
| `st5`    | `CoreNozzle`           | 5        | Core nozzle throat                 |
| `st6`    | `CoreNozzleExit`       | 6        | Core nozzle exit / downstream      |
| `st7`    | `FanNozzle`            | 7        | Fan nozzle throat                  |
| `st8`    | `FanNozzleExit`        | 8        | Fan nozzle exit / downstream       |
| `st9`    | `OfftakeDisch`         | 9        | Offtake air discharge              |

## Station-shortcut access

In addition to `eng.stN` (which returns the full `FlowStation`), every
thermodynamic quantity can be accessed via a `quantity + station_number`
shortcut:

```julia
eng.Tt4    # ≡ eng.st4.Tt  (total temperature at combustor exit)
eng.pt25c  # ≡ eng.st25c.pt (total pressure after inter-cooler)
eng.A9     # ≡ eng.st9.A   (area at offtake discharge)
eng.mdot2  # ≡ eng.st2.mdot
```

Shortcut reads and writes resolve to the same machine code as the
equivalent explicit nested access: the `Val{symbol}` dispatch
constant-folds at compile time.

## Constructors

```julia
EngineState{T}()   # all stations + ambient + design zeroed, numeric type T
EngineState()      # Float64 default
```
"""
mutable struct EngineState{T<:AbstractFloat}
    # -----------------------------------------------------------------------
    # Flow stations — one FlowStation per EngineStation member
    # Named by TASOPT station number with "st" prefix (alphabetical order
    # matches the EngineStation enumeration order for readability).
    # -----------------------------------------------------------------------
    st0   ::FlowStation{T}   # Freestream           — station 0
    st18  ::FlowStation{T}   # FanFaceOuter         — station 18
    st19  ::FlowStation{T}   # FanFaceLPC           — station 19
    st19c ::FlowStation{T}   # PreCoolerOut         — station 19c
    st2   ::FlowStation{T}   # FanFaceFan           — station 2
    st21  ::FlowStation{T}   # FanExit              — station 21
    st25  ::FlowStation{T}   # LPCExit              — station 25
    st25c ::FlowStation{T}   # InterCoolerOut       — station 25c
    st3   ::FlowStation{T}   # HPCExit              — station 3
    st4   ::FlowStation{T}   # CombustorExit        — station 4
    st4a  ::FlowStation{T}   # CoolMixInlet         — station 4a
    st41  ::FlowStation{T}   # TurbineInlet         — station 41
    st45  ::FlowStation{T}   # HPTExit              — station 45
    st49  ::FlowStation{T}   # LPTExit              — station 49
    st49c ::FlowStation{T}   # RegenCoolerOut       — station 49c
    st5   ::FlowStation{T}   # CoreNozzle           — station 5
    st6   ::FlowStation{T}   # CoreNozzleExit       — station 6
    st7   ::FlowStation{T}   # FanNozzle            — station 7
    st8   ::FlowStation{T}   # FanNozzleExit        — station 8
    st9   ::FlowStation{T}   # OfftakeDisch         — station 9

    # -----------------------------------------------------------------------
    # Ambient / flight-condition scalars
    # -----------------------------------------------------------------------
    M0 ::T   # freestream Mach number              [—]
    T0 ::T   # freestream static temperature       [K]
    p0 ::T   # freestream static pressure          [Pa]
    a0 ::T   # freestream speed of sound           [m/s]

    # -----------------------------------------------------------------------
    # Frozen design-point state (set by tfsize!, read by tfoper!)
    # -----------------------------------------------------------------------
    design ::DesignState{T}
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    EngineState{T}() where {T<:AbstractFloat}

Return a fully zero-initialised `EngineState` with numeric type `T`.
Every `FlowStation` is zeroed; all ambient scalars are `zero(T)`; the
embedded `DesignState` is zeroed.
"""
function EngineState{T}() where {T<:AbstractFloat}
    z  = zero(T)
    fs = FlowStation{T}
    EngineState{T}(
        fs(), fs(), fs(), fs(), fs(),   # st0, st18, st19, st19c, st2
        fs(), fs(), fs(), fs(), fs(),   # st21, st25, st25c, st3, st4
        fs(), fs(), fs(), fs(), fs(),   # st4a, st41, st45, st49, st49c
        fs(), fs(), fs(), fs(), fs(),   # st5, st6, st7, st8, st9
        z, z, z, z,                     # M0, T0, p0, a0
        DesignState{T}(),               # design
    )
end

"""
    EngineState()

Return a zero-initialised `EngineState{Float64}`.
"""
EngineState() = EngineState{Float64}()

# ---------------------------------------------------------------------------
# Station-shortcut property access — eng.Tt4 → eng.st4.Tt
# ---------------------------------------------------------------------------

"""
    _ENGINE_OWN_FIELDS

Tuple of every structural field name of `EngineState`.  Used by the
property-access overloads to distinguish direct field access from
station-shortcut resolution.
"""
const _ENGINE_OWN_FIELDS = (
    :st0, :st18, :st19, :st19c, :st2, :st21, :st25, :st25c,
    :st3, :st4, :st4a, :st41, :st45, :st49, :st49c,
    :st5, :st6, :st7, :st8, :st9,
    :M0, :T0, :p0, :a0, :design,
)

# Every FlowStation field name that is forwarded (all GasState fields + A + mdot).
const _FS_FORWARDED_FIELDS = (_GAS_FIELDS..., :A, :mdot)

# All shortcut symbols — used by propertynames() for tab-completion.
const _STATION_SHORTCUT_NAMES = let names = Symbol[]
    for stfld in _ENGINE_OWN_FIELDS
        s = String(stfld)
        startswith(s, "st") || continue
        suffix = s[3:end]   # "st4" → "4", "st19c" → "19c"
        for qty in _FS_FORWARDED_FIELDS
            push!(names, Symbol(qty, suffix))
        end
    end
    Tuple(names)
end

# ---------------------------------------------------------------------------
# Per-symbol generated method specialisations
#
# Using pure Val{symbol} dispatch eliminates branches in getproperty /
# setproperty!, giving the type inference engine a single concrete return
# type per call site.  The pattern is:
#
#   Base.getproperty(eng, :Tt4)
#     → _eng_getprop(eng, Val(:Tt4))     [one Val method per symbol]
#     → getproperty(getfield(eng, :st4), :Tt)
#     → Float64                           [fully inferred, no union]
#
# We generate methods for:
#   1. Every structural own field   — _eng_getprop returns getfield(eng, fld)
#   2. Every station shortcut       — _eng_getprop returns the nested quantity
# Setproperty follows the same pattern.
# ---------------------------------------------------------------------------

# Fallback — fires for any symbol not matched by the generated specialisations.
@inline _eng_getprop(::EngineState, ::Val{name}) where name =
    error("EngineState has no property $name")
@inline _eng_setprop!(::EngineState, ::Val{name}, _) where name =
    error("EngineState has no property $name")

let
    # ---- own structural fields ----
    for fld in _ENGINE_OWN_FIELDS
        @eval @inline _eng_getprop(eng::EngineState, ::Val{$(QuoteNode(fld))}) =
            getfield(eng, $(QuoteNode(fld)))
        @eval @inline _eng_setprop!(eng::EngineState, ::Val{$(QuoteNode(fld))}, v) =
            setfield!(eng, $(QuoteNode(fld)), v)
    end

    # ---- station-shortcut symbols ----
    for stfld in _ENGINE_OWN_FIELDS
        s = String(stfld)
        startswith(s, "st") || continue
        suffix = s[3:end]
        for qty in _FS_FORWARDED_FIELDS
            sc = Symbol(qty, suffix)
            @eval @inline _eng_getprop(eng::EngineState, ::Val{$(QuoteNode(sc))}) =
                getproperty(getfield(eng, $(QuoteNode(stfld))), $(QuoteNode(qty)))
            @eval @inline _eng_setprop!(eng::EngineState, ::Val{$(QuoteNode(sc))}, v) =
                setproperty!(getfield(eng, $(QuoteNode(stfld))), $(QuoteNode(qty)), v)
        end
    end
end

# ---------------------------------------------------------------------------
# Base overloads — thin wrappers that delegate to the Val-dispatched helpers
# ---------------------------------------------------------------------------

"""
    Base.getproperty(eng::EngineState, name::Symbol)

Return the named property of `eng`.  Direct struct fields (e.g. `eng.st4`,
`eng.M0`) and station-shortcut symbols (e.g. `eng.Tt4`, `eng.pt25c`,
`eng.A9`) are all resolved via `_eng_getprop(eng, Val(name))`.  Each
`Val{symbol}` specialisation has a single concrete return type, so the
compiler can fully infer the result of any literal-symbol access.
"""
@inline Base.getproperty(eng::EngineState, name::Symbol) =
    _eng_getprop(eng, Val(name))

"""
    Base.setproperty!(eng::EngineState, name::Symbol, v)

Set the named property of `eng`.  Direct struct fields and station-shortcut
symbols (e.g. `eng.Tt4 = 1600.0`) are all routed via
`_eng_setprop!(eng, Val(name), v)`.
"""
@inline Base.setproperty!(eng::EngineState, name::Symbol, v) =
    _eng_setprop!(eng, Val(name), v)

"""
    Base.propertynames(eng::EngineState, private::Bool=false)

Return all accessible property names: structural fields followed by every
station-shortcut symbol (e.g. `:Tt4`, `:pt25c`, `:A9`).
"""
Base.propertynames(::EngineState, private::Bool=false) =
    (_ENGINE_OWN_FIELDS..., _STATION_SHORTCUT_NAMES...)
