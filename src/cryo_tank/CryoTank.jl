"""
`CryoTank` is a module that contains the structural, thermal and energy models of a cryogenic tank.
"""
module CryoTank

using ..engine
using ..atmosphere
using ..structures
using ..TASOPT: fuselage_tank, aircraft
using ..materials
using Roots
using NLopt
using StaticArrays

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, compare_strings, FuelType, newton_solve, ode_solve

export SaturatedMixture, SaturatedPhaseProps, tanksize!,
       InnerTankGeometry, OuterTankGeometry

include(__TASOPTindices__)
include("../utils/constants.jl")

include("tankWmech.jl")
include("tankWthermal.jl")
include("tanksize.jl")

include("fuel_thermo.jl")
include("mixture.jl")
include("pressure.jl")
include("tanktools.jl")

end