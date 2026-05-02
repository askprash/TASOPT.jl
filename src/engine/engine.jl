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

# ── Public API: typed state and engine entry points ──────────────────────────
# These are the intended extension surface for aircraft sizing and analysis.
# Component equation kernels (compressor_pratd, turbine_efficiency, etc.) are
# exported so callers can assemble custom Newton systems; low-level gas functions
# and turbomachinery map helpers are internal — access them as engine.gassum etc.
export Engine

# State types and configuration
export CalcMode, CoolingOpt, EngineStation
export DesignState
export GasState, AIR_ALPHA
export FlowStation
export EngineState

# Component types and their Newton-kernel entry points
export Inlet, inlet_diffuser!, inlet_bli_mixing!
export Combustor, combustor_exit!, combustor_burnd
export Compressor, compressor_efficiency, compressor_exit!, compressor_Nb_residual, compressor_pratd
export Shaft, hp_shaft_work, lp_shaft_work, hp_shaft_workd, lp_shaft_workd, shaft_speed_residual
export Splitter, bypass_ratio
export Nozzle, nozzle_exit, nozzle_massflow_residual, nozzle_gross_thrust
export TurbineMap, Turbine, turbine_efficiency, turbine_delhd, turbine_exit!, turbine_mb_residual
export dump_stations

# Turbofan calculation entry points and results
export tfwrap!, tfcalc!, mcool
export tfweightwrap!, tfweight, tfsize!, tfoper!
export SizingResult

# Turbofan design-point and sweep runners
export run_engine_design_point
export SweepResult, run_engine_sweep, write_sweep_csv
export write_sweep_toml, ENGINE_BASELINE_PATH, regenerate_engine_baseline

# Ducted-fan entry points and results
export ductedfanoper!, ductedfansize!, ductedfancalc!, ductedfanweight, fuel_cell_with_ducted_fan_weight!
export DuctedFanState, engine_state_to_ducted_fan_state!
export run_ducted_fan_design_point, run_ducted_fan_sweep
export write_ducted_fan_sweep_toml, DUCTED_FAN_BASELINE_PATH, regenerate_ducted_fan_baseline

# Gas thermodynamics utilities used outside the engine module
export fuelLHV, gasPr

# Heat-exchanger public API
export hxdesign!, radiator_design!, hxweight, resetHXs, HXOffDesign!, RadiatorOffDesign!, check_HX_overwriting

# Fuel-cell / ducted-fan coupling
export calculate_fuel_cell_with_ducted_fan!, ductedfanweight!

export check_engine_convergence_failure

# ── Internal symbols (not exported) ──────────────────────────────────────────
# Access via engine.<symbol> when needed from tests or scripts:
#   station_number, station_description
#   set_total_from_Tt!, set_static_from_M!, apply_pratio_from!, apply_delh_from!
#   Tmcalc, gas_tset, gaschem
#   gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd
#   gas_mach, gas_machd, gas_mass, gasfuel
#   ddct, ddat, gct, gat
#   Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, StructuralAlloy, unpack_ac, compare_strings
import ..TASOPT.atmosphere: atmos

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils/constants.jl"))
include("turbofan/engine_enums.jl")
include("turbofan/design_state.jl")
include("turbofan/gas_state.jl")
include("turbofan/flow_station.jl")
include("turbofan/engine_state.jl")
include("turbofan/inlet.jl")
include("gasfun.jl")
include("gascalc.jl")
include("turbofan/thermo_wrappers.jl")
# include("tfan.jl")
include("turbomachinery/tfmap.jl")
include("turbomachinery/maps.jl")
include("turbofan/turbine.jl")
include("turbofan/combustor.jl")
include("turbofan/shaft.jl")
include("turbofan/compressor.jl")
include("turbofan/splitter.jl")
include("turbofan/nozzle.jl")
include("turbofan/tfcool.jl")
include("turbofan/sizing_result.jl")
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
include("ducted_fan/ducted_fan_harness.jl")

#Fuel cell models
include("PEMfuelcell.jl")
include("fuel_cell/FC_objects.jl")
include("fuel_cell/fuel_cell_operations.jl")
include("fuel_cell/FC_ducted_fan_models.jl")
include("fuel_cell/FC_ducted_fan_weight.jl")

end