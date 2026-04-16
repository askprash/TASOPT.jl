"""
Typed containers for per-flight-point and per-mission engine state,
replacing the bare `pare[:,ip,im]` array dimensions incrementally.

The hierarchy mirrors the legacy array layout:

    ac.missions                          :: Vector{Mission{Float64}}
      └─ missions[im]                    :: Mission{T}
           └─ points[ip]                 :: MissionPoint{T}
                └─ engine                :: EngineState{T}
                └─ (future) aero         :: AeroState{T}  (replaces para[:,ip,im])

During the current adapter-walk wave, `pare` remains the source of truth.
`pare_to_engine_state!` mirrors engine state into each `MissionPoint.engine`
after every `tfwrap!` call.  `pare` indices are retired incrementally as
downstream consumers migrate to typed access.
"""

# ---------------------------------------------------------------------------
# MissionPoint
# ---------------------------------------------------------------------------

"""
    MissionPoint{T<:AbstractFloat}

Per-flight-profile-point container.  One `MissionPoint` exists for each of
the `npoints` (typically `iptotal = 17`) operating points in a mission.

## Fields

| Field    | Type            | Description                                      |
|:---------|:----------------|:-------------------------------------------------|
| `engine` | `EngineState{T}`| Complete engine aerothermodynamic state          |

Future fields (`aero::AeroState{T}`, etc.) will replace the corresponding
`para[:,ip,im]` slice as the adapter-walk progresses.

## Constructors

```julia
MissionPoint{T}()   # zeroed engine state, numeric type T
MissionPoint()      # Float64 default
```
"""
mutable struct MissionPoint{T<:AbstractFloat}
    engine ::EngineState{T}
    # Future: aero ::AeroState{T}   (replaces para[:,ip,im] in a later wave)
end

MissionPoint{T}() where {T} = MissionPoint{T}(EngineState{T}())
MissionPoint()              = MissionPoint{Float64}()

# ---------------------------------------------------------------------------
# Mission
# ---------------------------------------------------------------------------

"""
    Mission{T<:AbstractFloat}

One complete flight mission: a `Vector` of `MissionPoint{T}`, one per
operating point.  The vector length equals `npoints` (typically `iptotal`).

## Fields

| Field    | Type                       | Description                             |
|:---------|:---------------------------|:----------------------------------------|
| `points` | `Vector{MissionPoint{T}}`  | Per-point state; length == npoints      |

## Constructors

```julia
Mission{T}(npoints::Int)   # npoints zero-initialised MissionPoints, type T
Mission(npoints::Int)      # Float64 default
Mission{T}()               # empty (length 0); populated by read_input
Mission()                  # Float64 empty
```
"""
mutable struct Mission{T<:AbstractFloat}
    points ::Vector{MissionPoint{T}}   # length == npoints (typically iptotal = 17)
end

Mission{T}(npoints::Int) where {T<:AbstractFloat} = Mission{T}([MissionPoint{T}() for _ in 1:npoints])
Mission(npoints::Int)                            = Mission{Float64}(npoints)
Mission{T}() where {T<:AbstractFloat}            = Mission{T}(0)   # empty; populated by read_input
Mission()                                        = Mission{Float64}()
