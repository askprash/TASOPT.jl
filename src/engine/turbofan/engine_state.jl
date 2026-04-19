"""
Engine-state container: all named flow stations, ambient flight condition,
and a composed design-state for a TASOPT turbofan.
"""

using Printf

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
    M0    ::T   # freestream Mach number              [—]
    T0    ::T   # freestream static temperature       [K]
    p0    ::T   # freestream static pressure          [Pa]
    a0    ::T   # freestream speed of sound           [m/s]
    rho0  ::T   # freestream density                  [kg/m³]
    mu0   ::T   # freestream dynamic viscosity        [Pa·s]
    Tfuel      ::T   # fuel temperature at combustor inlet [K] → ieTfuel
    Tfuel_tank ::T   # fuel temperature in tank (source of truth) [K] → ieTft
    hvap  ::T   # fuel latent heat of vaporisation (initial, const) [J/kg] → iehvap

    # -----------------------------------------------------------------------
    # Heat-exchanger delta outputs (tasopt-7vz / tasopt-j9l.41.2 / tasopt-w82)
    # Written to EngineState by HXOffDesign!/resetHXs.
    # Read by tfcalc! via the eng_hx parameter (tasopt-dti).
    # Legacy pare slots (iePreCDeltah etc.) removed in tasopt-w82.
    # -----------------------------------------------------------------------
    PreCDeltah    ::T   # pre-cooler enthalpy delta     [J/kg]
    PreCDeltap    ::T   # pre-cooler pressure delta     [Pa]
    InterCDeltah  ::T   # inter-cooler enthalpy delta   [J/kg]
    InterCDeltap  ::T   # inter-cooler pressure delta   [Pa]
    RegenDeltah   ::T   # regen-cooler enthalpy delta   [J/kg]
    RegenDeltap   ::T   # regen-cooler pressure delta   [Pa]
    TurbCDeltah   ::T   # turbine-cooler enthalpy delta [J/kg]
    TurbCDeltap   ::T   # turbine-cooler pressure delta [Pa]
    RadiatorDeltah ::T  # radiator enthalpy delta       [J/kg]
    RadiatorDeltap ::T  # radiator pressure delta       [Pa]
    HXrecircP     ::T   # HX recirculation pump power   [W]
    hvapcombustor ::T   # effective hvap in combustor   [J/kg]

    # -----------------------------------------------------------------------
    # Engine-level performance outputs
    # Written by tfcalc! EXIT blocks and flushed to pare via
    # engine_state_to_pare!.  Read back from pare via pare_to_engine_state!.
    # -----------------------------------------------------------------------
    TSFC  ::T   # thrust-specific fuel consumption [kg/N/s]
    Fe    ::T   # net thrust per engine            [N]
    Fsp   ::T   # specific thrust                  [N/(kg/s)]
    BPR   ::T   # bypass ratio                     [—]
    mfuel ::T   # total fuel mass flow (all engines) [kg/s]

    # -----------------------------------------------------------------------
    # Ducted-fan performance outputs (tasopt-j9l.45.3)
    # Written by ductedfancalc! EXIT block; backed by pare[iePfan], pare[ieTSEC],
    # pare[iemfan].  Zero for turbofan missions.
    # -----------------------------------------------------------------------
    Pfan  ::T   # fan shaft power                  [W]
    TSEC  ::T   # thrust-specific energy consumption [J/N]
    mfan  ::T   # fan total mass flow              [kg/s]

    # -----------------------------------------------------------------------
    # Compressor map operating points (tasopt-drd)
    # Written by tfcalc! EXIT block; read back via pare_to_engine_state!.
    # Synced to pare before off-design calls via engine_state_to_pare!.
    # Back-propagated from previous mission point as Newton initial guesses
    # (see _mission_iteration! descent initialisation).
    # -----------------------------------------------------------------------
    mbf  ::T   # fan corrected mass flow   [kg/s corrected]
    mblc ::T   # LPC corrected mass flow   [kg/s corrected]
    mbhc ::T   # HPC corrected mass flow   [kg/s corrected]
    pif  ::T   # fan pressure ratio        [—]
    pilc ::T   # LPC pressure ratio        [—]
    pihc ::T   # HPC pressure ratio        [—]

    # -----------------------------------------------------------------------
    # Per-point nozzle area schedule factors (tasopt-dw7)
    # A5fac: core nozzle throat area relative to design value  [—]
    # A7fac: fan  nozzle throat area relative to design value  [—]
    # Set by read_input.jl for each flight point; synced from
    # pare[ieA5fac] / pare[ieA7fac] via pare_to_engine_state!.
    # -----------------------------------------------------------------------
    A5fac   ::T   # core nozzle area factor             [—]
    A7fac   ::T   # fan  nozzle area factor             [—]
    Pfanmax ::T   # maximum fan shaft power (per-point input) [W]

    # -----------------------------------------------------------------------
    # Component adiabatic efficiencies (tasopt-j9l.63.1)
    # Written by tfcalc! EXIT blocks; backed by pare[ieetaf..ieetalt].
    # -----------------------------------------------------------------------
    etaf  ::T   # fan adiabatic efficiency           [—]
    etalc ::T   # LPC adiabatic efficiency           [—]
    etahc ::T   # HPC adiabatic efficiency           [—]
    etaht ::T   # HPT adiabatic efficiency           [—]
    etalt ::T   # LPT adiabatic efficiency           [—]

    # -----------------------------------------------------------------------
    # Overall propulsion efficiencies (tasopt-j9l.63.2)
    # Derived quantities — no pare backing; set to zero by pare_to_engine_state!.
    # Written by tfcalc! EXIT blocks only; engine_state_to_pare! is a no-op.
    #
    # Definitions (standard propulsion, Cumpsty Ch.1, Hill & Peterson Ch.5):
    #   eta_overall  = Fe * u0 / (mdotf_per_engine * hfuel)
    #   eta_prop     = 2 / (2 + Fsp)      [mixed-jet momentum: u_jet = u0*(1+Fsp)]
    #   eta_thermal  = eta_overall / eta_prop
    #
    # All three are zero at ground-idle (Fe ≤ 0 or Fsp ≤ 0).
    # -----------------------------------------------------------------------
    eta_thermal ::T   # thermal efficiency    [—]
    eta_prop    ::T   # propulsive efficiency [—]
    eta_overall ::T   # overall efficiency    [—]

    # -----------------------------------------------------------------------
    # Per-point operating outputs (tasopt-j9l.45.1)
    # Written by tfcalc! common EXIT block; dual-write alongside pare[ie*].
    # Backed by pare[ieConvFail], pare[iehfuel..ieKinl], pare[ieNf..ieNbhc],
    # pare[ieepf..ieeplt].
    # -----------------------------------------------------------------------
    ConvFail ::T   # convergence failure flag: 0.0 = converged, 1.0 = failed [—]
    hfuel    ::T   # fuel specific enthalpy                                    [J/kg]
    ff       ::T   # fuel-to-air ratio                                         [—]
    mofft    ::T   # mass offtake per engine                                   [kg/s]
    Pofft    ::T   # power offtake per engine                                  [W]
    Phiinl   ::T   # BLI ingested momentum-excess flux per engine              [W]
    Kinl     ::T   # BLI ingested kinetic-energy flux per engine               [W]

    # -----------------------------------------------------------------------
    # Spool speeds (tasopt-j9l.45.1)
    # Written by tfcalc! common EXIT block.
    # Backed by pare[ieNf], pare[ieN1], pare[ieN2], pare[ieNbf..ieNbhc].
    # -----------------------------------------------------------------------
    Nf   ::T   # fan physical spool speed (normalised by design)           [—]
    N1   ::T   # LP spool physical speed (normalised by design)            [—]
    N2   ::T   # HP spool physical speed (normalised by design)            [—]
    Nbf  ::T   # fan corrected spool speed (normalised)                    [—]
    Nblc ::T   # LPC corrected spool speed (normalised)                    [—]
    Nbhc ::T   # HPC corrected spool speed (normalised)                    [—]

    # -----------------------------------------------------------------------
    # Component polytropic loss fractions (tasopt-j9l.45.1)
    # Returned by tfsize!/tfoper!; written by tfcalc! common EXIT block.
    # Backed by pare[ieepf..ieeplt].
    # -----------------------------------------------------------------------
    epf  ::T   # fan polytropic loss fraction                              [—]
    eplc ::T   # LPC polytropic loss fraction                              [—]
    ephc ::T   # HPC polytropic loss fraction                              [—]
    epht ::T   # HPT polytropic loss fraction                              [—]
    eplt ::T   # LPT polytropic loss fraction                              [—]

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
        z, z, z, z, z, z, z, z, z,      # M0, T0, p0, a0, rho0, mu0, Tfuel, Tfuel_tank, hvap
        z, z, z, z, z, z, z, z,        # HX Δh/Δp: PreC, InterC, Regen, TurbC
        z, z, z, z,                     # HX Δh/Δp: Radiator, recircP, hvapcombustor
        z, z, z, z, z,                  # TSFC, Fe, Fsp, BPR, mfuel
        z, z, z,                        # Pfan, TSEC, mfan
        z, z, z, z, z, z,              # mbf, mblc, mbhc, pif, pilc, pihc
        z, z, z,                        # A5fac, A7fac, Pfanmax
        z, z, z, z, z,                  # etaf, etalc, etahc, etaht, etalt
        z, z, z,                        # eta_thermal, eta_prop, eta_overall
        z, z, z, z, z, z, z,           # ConvFail, hfuel, ff, mofft, Pofft, Phiinl, Kinl
        z, z, z, z, z, z,              # Nf, N1, N2, Nbf, Nblc, Nbhc
        z, z, z, z, z,                  # epf, eplc, ephc, epht, eplt
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
    :M0, :T0, :p0, :a0, :rho0, :mu0, :Tfuel, :Tfuel_tank, :hvap,
    :PreCDeltah, :PreCDeltap, :InterCDeltah, :InterCDeltap,
    :RegenDeltah, :RegenDeltap, :TurbCDeltah, :TurbCDeltap,
    :RadiatorDeltah, :RadiatorDeltap, :HXrecircP, :hvapcombustor,
    :TSFC, :Fe, :Fsp, :BPR, :mfuel,
    :Pfan, :TSEC, :mfan,
    :mbf, :mblc, :mbhc, :pif, :pilc, :pihc,
    :A5fac, :A7fac, :Pfanmax,
    :etaf, :etalc, :etahc, :etaht, :etalt,
    :eta_thermal, :eta_prop, :eta_overall,
    :ConvFail, :hfuel, :ff, :mofft, :Pofft, :Phiinl, :Kinl,
    :Nf, :N1, :N2, :Nbf, :Nblc, :Nbhc,
    :epf, :eplc, :ephc, :epht, :eplt,
    :design,
)

# FlowStation fields forwarded as EngineState station-shortcut symbols.
# NOTE: :st (total entropy complement) is deliberately excluded — the symbol
# prefix "st" is reserved for station fields (st0, st18, st2, …), so generating
# shortcuts like Symbol(:st, "4") == :st4 would silently shadow the structural
# station field.  Access station entropy via the explicit path: eng.st4.st
const _FS_FORWARDED_FIELDS = (
    :Tt, :ht, :pt, :cpt, :Rt,           # total thermo  (:st excluded — see note above)
    :Ts, :ps, :hs, :ss, :cps, :Rs, :u,  # static thermo
    :alpha,                               # composition
    :A, :mdot,                           # station mechanics
)

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

    # ---- derived shortcuts (tasopt-j9l.63.1) ----
    # eng.mcore → eng.st2.mdot  (core mass flow; no separate field)
    @eval @inline _eng_getprop(eng::EngineState, ::Val{:mcore}) =
        getproperty(getfield(eng, :st2), :mdot)
    @eval @inline _eng_setprop!(eng::EngineState, ::Val{:mcore}, v) =
        setproperty!(getfield(eng, :st2), :mdot, v)
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

# ---------------------------------------------------------------------------
# Station dump utility
# ---------------------------------------------------------------------------

# Physical flow-path order: (TASOPT#, short name, EngineState field symbol)
const _STATION_DUMP_ORDER = (
    ("0",   "Freestream",    :st0),
    ("2",   "FanFaceFan",    :st2),
    ("18",  "FanFaceOuter",  :st18),
    ("19",  "FanFaceLPC",    :st19),
    ("19c", "PreCoolerOut",  :st19c),
    ("21",  "FanExit",       :st21),
    ("25",  "LPCExit",       :st25),
    ("25c", "InterCoolerOut",:st25c),
    ("3",   "HPCExit",       :st3),
    ("4",   "CombustorExit", :st4),
    ("4a",  "CoolMixInlet",  :st4a),
    ("41",  "TurbineInlet",  :st41),
    ("45",  "HPTExit",       :st45),
    ("49",  "LPTExit",       :st49),
    ("49c", "RegenCoolerOut",:st49c),
    ("5",   "CoreNozzle",    :st5),
    ("6",   "CoreNozzleExit",:st6),
    ("7",   "FanNozzle",     :st7),
    ("8",   "FanNozzleExit", :st8),
    ("9",   "OfftakeDisch",  :st9),
)

"""
    dump_stations([io::IO,] eng::EngineState)

Print a human-readable table of all twenty named flow stations in `eng`.
Rows are stations in physical flow order; columns are total state
(Tt, pt, ht), static state (Ts, ps, u), area (A), and mass flow (mdot).

# Arguments
- `io`: output target; defaults to `stdout`.
- `eng`: engine state to inspect.

# Example
```julia
dump_stations(eng)                    # prints to stdout
dump_stations(stderr, eng)            # prints to stderr
open("stations.txt", "w") do f
    dump_stations(f, eng)
end
```
"""
function dump_stations(io::IO, eng::EngineState)
    @printf(io, "Engine stations — M0=%.4f  T0=%.2fK  p0=%.2fPa  a0=%.2fm/s\n",
            eng.M0, eng.T0, eng.p0, eng.a0)
    println(io)
    @printf(io, "%-4s  %-14s  %9s  %12s  %12s  %9s  %12s  %8s  %9s  %11s\n",
            "St#", "Name",
            "Tt[K]", "pt[Pa]", "ht[J/kg]",
            "Ts[K]", "ps[Pa]", "u[m/s]",
            "A[m²]", "mdot[kg/s]")
    println(io, "─"^108)
    for (stnum, stname, stfld) in _STATION_DUMP_ORDER
        fs = getfield(eng, stfld)
        @printf(io, "%-4s  %-14s  %9.2f  %12.5g  %12.5g  %9.2f  %12.5g  %8.2f  %9.5f  %11.5f\n",
                stnum, stname,
                fs.Tt, fs.pt, fs.ht,
                fs.Ts, fs.ps, fs.u,
                fs.A, fs.mdot)
    end
end

"""
    dump_stations(eng::EngineState)

Print station table to `stdout`.  See `dump_stations(io, eng)` for full docs.
"""
dump_stations(eng::EngineState) = dump_stations(stdout, eng)
