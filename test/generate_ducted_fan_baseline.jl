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
parg = ac.parg

enginecalc! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
engineweight! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
enginemodel = TASOPT.engine.FuelCellDuctedFan(
    "fuel_cell_with_ducted_fan", enginecalc!, "nasa", engineweight!, false)

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
ac.para[iaROCdes, ipclimb1:ipclimbn, :] .= 500 * TASOPT.ft_to_m / 60

engine = TASOPT.engine.Engine(enginemodel, fcdata, Vector{TASOPT.engine.HeatExchanger}())
ac.engine = engine

# Seed per-point typed state: Pfanmax and radiator fields.
for im in 1:length(ac.missions), ip in 1:length(ac.missions[im].points)
    eng_ip = ac.missions[im].points[ip].engine
    eng_ip.Pfanmax     = 10e6
    eng_ip.RadCoolantT = engine.data.FC_temperature[ip, im]
    eng_ip.RadCoolantP = engine.data.FC_pressure[ip, im]
    eng_ip.Qradiator   = engine.data.FC_heat[ip, im]
end

# Design-point engine inputs at ipcruise1 written directly to typed state.
let eng = ac.missions[1].points[ipcruise1].engine
    eng.Fe           = 16981.808185580507
    eng.Fsp          = 0.5268888878557653
    eng.pif          = 1.685
    eng.design.pid   = 0.998
    eng.design.pifn  = 0.98
    eng.design.epolf = 0.8948
    eng.design.pifK  = 1.685
    eng.design.epfK  = -0.077
    eng.design.M2    = 0.6
end
parg[igneng]  = 2
ac.wing.layout.S = 81.25043040696103

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
