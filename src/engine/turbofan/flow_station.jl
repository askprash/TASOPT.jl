"""
Flow-station container: gas state plus cross-sectional area and mass flow.
"""

using StaticArrays

# ---------------------------------------------------------------------------
# GasState fields that FlowStation forwards to its embedded gas::GasState
# ---------------------------------------------------------------------------

# Tuple of every field name that lives in GasState.
# Declared as a constant so the compiler can specialise the `in` check
# to a constant fold at call sites.
const _GAS_FIELDS = (
    :Tt, :ht, :pt, :cpt, :Rt, :st,      # total thermodynamic state
    :Ts, :ps, :hs, :ss, :cps, :Rs, :u,  # static thermodynamic state
    :alpha,                               # species composition
)

# ---------------------------------------------------------------------------
# FlowStation
# ---------------------------------------------------------------------------

"""
    FlowStation{T<:AbstractFloat}

Mutable container for a single aerothermodynamic flow station.  Combines the
complete thermodynamic state (via an embedded `GasState{T}`) with the two
mechanical quantities that vary cross-section to cross-section in a duct:

| Field  | Unit  | Description                                          |
|:-------|:------|:-----------------------------------------------------|
| `gas`  | —     | Embedded `GasState{T}` (total + static thermo state) |
| `A`    | m²    | Cross-sectional flow area                            |
| `mdot` | kg/s  | Mass flow rate through the station                   |

## Ergonomic property forwarding

All thirteen `GasState` field names are forwarded transparently so callers
can write, e.g., `station.Tt` instead of `station.gas.Tt`.  `getproperty`
and `setproperty!` both honour the forwarding; the compiler constant-folds
the symbol dispatch so the generated machine code is identical to a direct
nested access.

## Constructors

Three constructors mirror the `GasState` constructors, with optional `A`
(area) and `mdot` (mass flow) keyword arguments that default to `zero(T)`:

```julia
FlowStation{T}(; A=zero(T), mdot=zero(T))         # all-zero gas + area + mdot
FlowStation(; A=0.0, mdot=0.0)                     # Float64 default
FlowStation{T}(Tt, ht, pt, cpt, Rt, alpha;         # explicit total state
               A=zero(T), mdot=zero(T))
```
"""
mutable struct FlowStation{T<:AbstractFloat}
    gas  ::GasState{T}   # complete thermodynamic state
    A    ::T             # cross-sectional area                  [m²]
    mdot ::T             # mass flow rate                        [kg/s]
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    FlowStation{T}(; A=zero(T), mdot=zero(T)) where {T<:AbstractFloat}

Return a zero-initialised `FlowStation` with numeric type `T`.
The embedded `GasState` is fully zeroed; `A` and `mdot` default to zero.
"""
function FlowStation{T}(; A::T=zero(T), mdot::T=zero(T)) where {T<:AbstractFloat}
    FlowStation{T}(GasState{T}(), A, mdot)
end

"""
    FlowStation(; A=0.0, mdot=0.0)

Return a zero-initialised `FlowStation{Float64}`.
"""
FlowStation(; A::Float64=0.0, mdot::Float64=0.0) = FlowStation{Float64}(; A, mdot)

"""
    FlowStation{T}(Tt, ht, pt, cpt, Rt, alpha; A=zero(T), mdot=zero(T))

Construct a `FlowStation` with explicit total-state fields and species
composition; static-state fields and momentum quantities default to zero.

Arguments mirror `GasState{T}(Tt, ht, pt, cpt, Rt, alpha)`.
"""
function FlowStation{T}(
        Tt::T, ht::T, pt::T, cpt::T, Rt::T, alpha::SVector{5,T};
        A::T=zero(T), mdot::T=zero(T),
) where {T<:AbstractFloat}
    FlowStation{T}(GasState{T}(Tt, ht, pt, cpt, Rt, alpha), A, mdot)
end

# ---------------------------------------------------------------------------
# Property forwarding — GasState fields accessible directly on FlowStation
# ---------------------------------------------------------------------------

"""
    Base.getproperty(fs::FlowStation, name::Symbol)

Forward any `GasState` field name to `fs.gas`; return `FlowStation` fields
(`gas`, `A`, `mdot`) directly.  The `name in _GAS_FIELDS` check constant-
folds at compile time, so no runtime dispatch overhead is incurred.
"""
@inline function Base.getproperty(fs::FlowStation, name::Symbol)
    name in _GAS_FIELDS && return getproperty(getfield(fs, :gas), name)
    return getfield(fs, name)
end

"""
    Base.setproperty!(fs::FlowStation, name::Symbol, val)

Forward assignment to any `GasState` field name to `fs.gas`; set
`FlowStation`-own fields (`A`, `mdot`) directly.
Assigning to `gas` replaces the embedded `GasState` as a whole.
"""
@inline function Base.setproperty!(fs::FlowStation, name::Symbol, val)
    if name in _GAS_FIELDS
        setproperty!(getfield(fs, :gas), name, val)
    else
        setfield!(fs, name, val)
    end
end

"""
    Base.propertynames(fs::FlowStation, private::Bool=false)

Return the union of `FlowStation` own fields and all forwarded `GasState`
field names.  Enables tab-completion and `propertynames` introspection.
"""
Base.propertynames(::FlowStation, private::Bool=false) =
    (:gas, :A, :mdot, _GAS_FIELDS...)
