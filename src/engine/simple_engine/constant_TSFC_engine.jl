""" 
    constant_TSFC_engine!(ac, d1, imission::Int64, ip::Int64, d2, d3)

A very simple model of an aircraft engine consuming some fuel. The fuel flow is determined by the thrust required 
by the aircraft. The thrust is calculated based on the climb rate of the aircraft. The fuel flow is calculated using 
the thrust specific fuel consumption (TSFC) of the engine, which is an input.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `d`s: dummy variables used for compatibility with other sizing functions

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function constant_TSFC_engine!(ac, d1, imission::Int64, ip::Int64, d2, d3 = 0)
    #Unpack data storage arrays
    parg, _, para, _, _, _, wing, _, _, _, _ = unpack_ac(ac, imission, ip = ip)
    eng = ac.missions[imission].points[ip].engine
    TSFC = eng.TSFC #Extract TSFC at this mission point
    neng = parg[igneng]

    if (ip in range(ipstatic, ipclimbn)) #If aircraft is in climb, thrust is determined from the climb rate

        calculate_thrust_from_ROC!(ac, ip, imission)
        Fe = eng.Fe #Extract the required thrust

    else
        Fe = eng.Fe #In other mission points, thrust is already computed

    end

    mfuel_per_eng =  TSFC * Fe / gee #Calculate fuel flow per engine from set TSFC
    eng.mfuel = mfuel_per_eng * neng
end