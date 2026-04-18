"""
generate_ducted_fan_baseline.jl — Regenerate the ducted-fan sweep regression
baseline at test/fixtures/ducted_fan_sweep_baseline.toml.

Run from the repository root:
    julia --project=. test/generate_ducted_fan_baseline.jl

Only run this script when the ducted-fan physics or test aircraft parameters
have intentionally changed and the new values have been verified.  Commit
the updated fixture file with an explanatory git message so reviewers can
audit the numerical diff.
"""

using TASOPT, TOML
include(TASOPT.__TASOPTindices__)

# ---------------------------------------------------------------------------
# Build the test ducted-fan aircraft — same configuration as the
# "Ducted fan harness" testset in test/unit_test_ductedfan.jl.
# ---------------------------------------------------------------------------

println("Building test ducted-fan aircraft...")

ac = TASOPT.load_default_model()
pare = ac.pare
parg = ac.parg

enginecalc! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
engineweight! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
enginemodel = TASOPT.engine.FuelCellDuctedFan(
    "fuel_cell_with_ducted_fan", enginecalc!, "nasa", engineweight!, false)
pare[iePfanmax,:,:] .= 10e6

fcdata = TASOPT.engine.FuelCellDuctedFanData(2)
fcdata.type = "HT-PEMFC"
fcdata.current_density[iprotate,:] .= 1e4
fcdata.FC_temperature .= 453.15
fcdata.FC_pressure .= 3e5
fcdata.water_concentration_anode .= 0.1
fcdata.water_concentration_cathode .= 0.1
fcdata.λ_H2 .= 3.0
fcdata.λ_O2 .= 3.0
fcdata.thickness_membrane = 100e-6
fcdata.thickness_anode    = 250e-6
fcdata.thickness_cathode  = 250e-6
fcdata.design_voltage = 200.0
pare[ieRadiatorepsilon,:,:] .= 0.7
pare[ieRadiatorMp,:,:]      .= 0.12
pare[ieDi,:,:]              .= 0.4
ac.para[iaROCdes, ipclimb1:ipclimbn, :] .= 500 * TASOPT.ft_to_m / 60

engine = TASOPT.engine.Engine(enginemodel, fcdata, Vector{TASOPT.engine.HeatExchanger}())
ac.engine = engine
pare[ieRadiatorCoolantT,:,:] = engine.data.FC_temperature[:,:]
pare[ieRadiatorCoolantP,:,:] = engine.data.FC_pressure[:,:]
pare[ieRadiatorHeat,:,:]     = engine.data.FC_heat[:,:]

pare[ieFe,    ipcruise1, 1] = 16981.808185580507
parg[igneng]  = 2
ac.wing.layout.S = 81.25043040696103
pare[ieFsp,   ipcruise1, 1] = 0.5268888878557653
pare[iepif,   ipcruise1, 1] = 1.685
pare[iepid,   ipcruise1, 1] = 0.998
pare[iepifn,  ipcruise1, 1] = 0.98
pare[ieepolf, ipcruise1, 1] = 0.8948
pare[iepifK,  ipcruise1, 1] = 1.685
pare[ieepfK,  ipcruise1, 1] = -0.077
pare[ieM2,    ipcruise1, 1] = 0.6

ac.para[iaalt,  ipcruise1, 1] = 10668.0   # FL350 standard day
ac.para[iaMach, ipcruise1, 1] = 0.8

# ---------------------------------------------------------------------------
# Run design point and single-point sweep.
# ---------------------------------------------------------------------------

println("Running design point and sweep...")

df = TASOPT.engine.run_ducted_fan_design_point(ac)
println("  Fe = $(df.Fe) N")
println("  pif = $(df.pif)")
println("  etaf = $(df.etaf)")

states = TASOPT.engine.run_ducted_fan_sweep(ac; ip_range=ipcruise1:ipcruise1)

# ---------------------------------------------------------------------------
# Write baseline.
# ---------------------------------------------------------------------------

TASOPT.engine.regenerate_ducted_fan_baseline(ac, ipcruise1:ipcruise1)
