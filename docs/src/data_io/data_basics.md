# [Data structure basics](@id datastructs_basics) 

Performance and design data is held largely in the `par` arrays (a holdover from FORTRAN TASOPT) along with a growing body of `struct`s to represent cohesive components and systems. An `aircraft` `struct` wraps these arrays, `struct`s, and auxiliary information.

!!! compat "Future Changes"
    We don't like this hybrid approach either. It's the legacy of Fortran.

    In a future major revision, we aim to completely replace the `par` array system with the `struct`-oriented approach.

    **Engine quantities are already migrated.** After any `size_aircraft!` or `fly_mission!` call, engine outputs for each flight point are accessible via the typed `EngineState` API (see [Engine typed state](@ref engine_typed_state) below) instead of `pare[ie*]`.

## `par` arrays

Three arrays contain both prescribed inputs *and* computed outputs, some of which are multi-dimensional:

  - **`parg`**`::AbstractVector{Float64}`:  **geometric** quantities, weights, structural values, and other values inherent to an instantiated design (not operation).
  - **`parm`**`::AbstractArray{Float64, 2}`:  **mission**-prescribing parameters. The second dimension allows the specification of multiple mission profiles.
  - **`para`**`::AbstractArray{Float64, 3}`:  **aerodynamic** performance quantities. The second dimension captures the variation over a mission. The third dimension allows the specification of multiple mission profiles.

Data in the `par` arrays are accessed via Integer indices defined at `src/data_structs/index.inc`. These indices can be added to a namespace via `include(__TASOPTindices__)`:

```julia
using TASOPT
#using __TASOPTroot__, which fetches the src directory
include(joinpath(__TASOPTroot__, "data_structs/index.inc"))

#or more concisely
include(__TASOPTindices__)
```

The variable names of these indices indicate which `par` array they should access and hint at the quantity in question. For example, `iaMach` retrieves the flight Mach number via `para[iaMach, ip, im]`, where `ip` is the flight-point index and `im` is the mission index.

Note that for the multi-dimensional `par` arrays, indexing with a single Integer only retrieves the value for the first flight point of the first mission (namely, the design mission). Additional indexing is required to access data from different flight points or missions. Indices for specific flight points are defined in `index.inc` and should be used when indexing `para`, e.g., `ipstatic` for static ground condition or `ipcruise1` for the start of cruise.


```@example dataaccess
using TASOPT
include(__TASOPTindices__)
ac = load_default_model()
size_aircraft!(ac)

# Single element — first flight point of the design mission
println("Single element: ", size(ac.para[iaMach]))
println(ac.para[iaMach])

# Full slice over all flight points and missions
println("All flight points, all missions: ", size(ac.para[iaMach,:,:]))

# Slice at a specific flight point across all missions
println("Mach at cruise start: ", size(ac.para[iaMach,ipcruise1,:]))
println(ac.para[iaMach,ipcruise1,:])

```

## [Engine typed state](@id engine_typed_state)

Engine performance outputs are also accessible through the typed `EngineState` API, which is the preferred way to read engine results after a solve. Each `FlightPoint` in a mission carries an `engine` field of type `EngineState`:

```julia
using TASOPT
include(__TASOPTindices__)
ac = load_default_model()
size_aircraft!(ac)

# Read cruise engine state (mission 1, cruise-start flight point)
eng = ac.missions[1].points[ipcruise1].engine

# Thermodynamic station quantities (total pressure, total temperature)
OPR  = eng.pt3 / eng.pt2        # Overall Pressure Ratio  (HPC exit / fan face)
FPR  = eng.pt21 / eng.pt2       # Fan Pressure Ratio      (fan exit / fan face)
Tt4  = eng.Tt4                  # Combustor exit temperature [K]
Tt3  = eng.Tt3                  # HPC exit temperature [K]

# Cycle-level outputs
TSFC = eng.TSFC                 # Thrust-specific fuel consumption [kg/N/s]
BPR  = eng.BPR                  # Bypass ratio
pif  = eng.pif                  # Fan pressure ratio
pihc = eng.pihc                 # HPC pressure ratio

# Aggregate over all flight points in a mission
Tt3_max = maximum(p.engine.Tt3 for p in ac.missions[1].points)
```

Station shortcuts follow the pattern `eng.TtN` / `eng.ptN` for total temperature and pressure at station `N` (e.g. `eng.Tt3`, `eng.pt21`). The full station list and all available fields are documented in [`EngineState`](@ref).






## `aircraft` `struct`

An `aircraft` is composed of `par` array fields, title and description fields, and a `is_sized` flag to indicate its status. An optional `fuse_tank` field is present as a trial for future `struct`-based development. All fields are dot-accessible and array elements can be changed (e.g., `ac.parg[igS] = 20`), though the `struct` itself is not mutable.

Refer to the [`struct` reference page](@ref datastructs) for add'l details.
