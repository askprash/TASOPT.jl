"""
`engine` is a module that contains all low-fidelity calculations
required in the aircraft sizing. 
"""
module engine

using NLopt
using Roots
using NLsolve
using LinearAlgebra
using Random

export Engine

export CalcMode, CoolingOpt, EngineStation, station_number, station_description
export DesignState
export GasState
export FlowStation
export EngineState
export dump_stations
export pare_to_engine_state!, run_engine_design_point
export SweepResult, run_engine_sweep, write_sweep_csv
export write_sweep_toml, ENGINE_BASELINE_PATH, regenerate_engine_baseline
export set_total_from_Tt!, set_static_from_M!, apply_pratio_from!, apply_delh_from!
export tfwrap!, tfcalc!, mcool, Tmcalc, gas_tset, gaschem
export tfweightwrap!, tfweight, ddct, ddat, gct, gat, tfsize!, Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap, tfoper!
export ductedfanoper!, ductedfansize!, ductedfancalc!, ductedfanweight, fuel_cell_with_ducted_fan_weight!

export gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd, gas_mach, gas_machd, gas_mass, gasfuel, fuelLHV, gasPr
export hxdesign!, radiator_design!, hxweight, resetHXs, HXOffDesign!, RadiatorOffDesign!, check_HX_overwriting
export calculate_fuel_cell_with_ducted_fan!, ductedfanweight!

export check_engine_convergence_failure

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, StructuralAlloy, unpack_ac, compare_strings
import ..TASOPT.atmosphere: atmos

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils/constants.jl"))
include("turbofan/engine_enums.jl")
include("turbofan/design_state.jl")
include("turbofan/gas_state.jl")
include("turbofan/flow_station.jl")
include("turbofan/engine_state.jl")
include("gasfun.jl")
include("gascalc.jl")
include("turbofan/thermo_wrappers.jl")
# include("tfan.jl")
include("turbomachinery/tfmap.jl")
include("turbomachinery/maps.jl")
include("turbofan/tfcool.jl")
include("turbofan/tfsize.jl")
include("thrust_from_ROC.jl")
include("gaussn.jl")
include("compare.jl")
include("turbofan/tfoper.jl")
include("turbofan/tfcalc.jl")
include("turbofan/tfweight.jl")
include("turbofan/tfwrap.jl")
include("turbofan/engine_harness.jl")
include("turbofan/tfweightwrap.jl")
include("hxfun.jl")
include("simple_engine/constant_TSFC_engine.jl")
include("simple_engine/fractional_engine_weight.jl")
include(joinpath(__TASOPTroot__,"data_structs/engine.jl"))
include("ducted_fan/ductedfancalc.jl")
include("ducted_fan/ductedfansize.jl")
include("ducted_fan/ductedfanoper.jl")
include("ducted_fan/ductedfanweight.jl")

#Fuel cell models
include("PEMfuelcell.jl")
include("fuel_cell/FC_objects.jl")
include("fuel_cell/fuel_cell_operations.jl")
include("fuel_cell/FC_ducted_fan_models.jl")
include("fuel_cell/FC_ducted_fan_weight.jl")

end