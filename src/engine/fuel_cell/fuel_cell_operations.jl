"""
    size_fuel_cell!(ac, ip, imission)

Sizes the fuel cell for the given mission point based on the design power, operating conditions, 
and fuel cell parameters. This function calls the `PEMsize` function to calculate the number of cells, 
cell area, and heat dissipation required to achieve the desired power output.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object containing engine and fuel cell data
    - `ip::Int64`: mission point index
    - `imission::Int64`: mission index

    **Output:**
    No direct outputs. The `ac.engine.data` object is modified with:
    - `fcdata.number_cells` → number of fuel cells required to meet power demand
    - `fcdata.area_cell` → area of each individual fuel cell
    - `fcdata.FC_heat[ip, imission]` → heat rejected by the fuel cell at the mission point
"""
function size_fuel_cell!(ac, ip, imission)
    fcdata = ac.engine.data
    
    Pdes = fcdata.design_power
    u_FC = engine.PEMFC_inputs()

    u_FC.j = fcdata.current_density[ip, imission]
    u_FC.T = fcdata.FC_temperature[ip, imission]	
    u_FC.p_A = fcdata.FC_pressure[ip, imission]
    u_FC.p_C = fcdata.FC_pressure[ip, imission]
    u_FC.x_H2O_A = fcdata.water_concentration_anode[ip, imission]
    u_FC.x_H2O_C = fcdata.water_concentration_cathode[ip, imission]
    u_FC.λ_H2 = fcdata.λ_H2[ip, imission]
    u_FC.λ_O2 = fcdata.λ_O2[ip, imission]
    u_FC.t_M = fcdata.thickness_membrane
    u_FC.t_A = fcdata.thickness_anode
    u_FC.t_C = fcdata.thickness_cathode
    u_FC.type = fcdata.type

    Vdes = fcdata.design_voltage
    #α_guess = fcdata.α_water[ip, imission]

    #n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_FC, α_guess)
    n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_FC)

    fcdata.FC_heat[ip, imission] = Q
    fcdata.number_cells = n_cells
    fcdata.area_cell = A_cell
end

"""
    operate_fuel_cell!(ac, ip, imission)

Operates the fuel cell at a given mission point by computing the power output, heat dissipation, 
hydrogen consumption, and current density based on the operating conditions and fuel cell parameters. 
This function uses the `PEMoper` function to simulate the fuel cell behavior.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object containing engine and fuel cell data
    - `ip::Int64`: mission point index
    - `imission::Int64`: mission index

    **Output:**
    No direct outputs. The following fields in the `ac` object are modified:
    - `engine.mfuel` → fuel mass flow rate [kg/s]
    - `fcdata.stack_voltage[ip, imission]` → fuel cell stack voltage [V]
    - `fcdata.FC_heat[ip, imission]` → heat rejected by the fuel cell [W]
    - `fcdata.current_density[ip, imission]` → current density [A/m²]
"""
function operate_fuel_cell!(ac, ip, imission)
    fcdata = ac.engine.data
    
    u_FC = engine.PEMFC_inputs()

    u_FC.T = fcdata.FC_temperature[ip, imission]	
    u_FC.p_A = fcdata.FC_pressure[ip, imission]
    u_FC.p_C = fcdata.FC_pressure[ip, imission]
    u_FC.x_H2O_A = fcdata.water_concentration_anode[ip, imission]
    u_FC.x_H2O_C = fcdata.water_concentration_cathode[ip, imission]
    u_FC.λ_H2 = fcdata.λ_H2[ip, imission]
    u_FC.λ_O2 = fcdata.λ_O2[ip, imission]
    u_FC.t_M = fcdata.thickness_membrane
    u_FC.t_A = fcdata.thickness_anode
    u_FC.t_C = fcdata.thickness_cathode
    u_FC.type = fcdata.type

    n_cells = fcdata.number_cells
    A_cell = fcdata.area_cell
    #j_guess = fcdata.current_density[ip, imission]
    #α_guess = fcdata.α_water[ip, imission]

    P = fcdata.FC_power[ip, imission]
    #mfuel, V_stack, Q, j, α = engine.PEMoper(P, n_cells, A_cell, u_FC)
    V_stack, Q = engine.PEMoper(P, n_cells, A_cell, u_FC)

    #Calculate fuel mass flow rate
    M_h2 = 2.016e-3 #kg/mol
    F = 96485.3329 # C/mol, Faraday contant
    I = P / V_stack #Current
    mfuel = I * M_h2 / (2*F) #Molar flux of hydrogen gas

    ac.missions[imission].points[ip].engine.mfuel = mfuel * ac.parg[igneng]
    fcdata.stack_voltage[ip, imission] = V_stack
    fcdata.FC_heat[ip, imission] = Q
    fcdata.current_density[ip, imission] = I/A_cell
    #fcdata.α_water[ip, imission] = α
  
end