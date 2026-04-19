"""
    calculate_fuel_cell_with_ducted_fan!(ac, case, imission, ip, initializes_engine, iterw = 0)

Calculates the performance of a fuel cell powering a ducted fan for the given mission point. 
This function handles both the design and off-design performance of the ducted fan and fuel cell system.

In the design case, the function sizes the ducted fan at the start of cruise and the fuel cell at takeoff conditions.
In the off-design case, it computes the power and thrust requirements for each mission segment.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object containing geometry, mission, and propulsion data
    - `case::String`: case identifier, e.g., "design" or "off_design"
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `initializes_engine::Bool`: flag to initialize the ducted fan and fuel cell model
        - `true`  → initialize variables for iteration
        - `false` → use current variables as initial guesses
    - `iterw::Int64`: sizing loop iteration (optional, defaults to 0)

    **Output:**
    No direct outputs. The `ac` object is modified with:
    - Updated fuel cell power (`fcdata.FC_power`)
    - Updated ducted fan performance parameters (`pare`)
    - Updated engine size and power requirements during the design phase
"""
function calculate_fuel_cell_with_ducted_fan!(ac, case, imission, ip, initializes_engine, iterw = 0)
    #Unpack aircraft data
    parg, _, para, pare, _, _, _, _, _, _, _, _ = unpack_ac(ac, imission)
    fcdata = ac.engine.data #Extract fuel cell data

    # Sync bare pare → typed EngineState for the current ip before ductedfancalc!
    # reads from it (tasopt-j9l.45.3). Callers that already sync (mission_iteration.jl)
    # pay a small redundant copy; callers that do not (standalone tests) need this.
    pare_to_engine_state!(ac.missions[imission].points[ip].engine,
                          view(ac.pare, :, ip, imission))
    eng_ip = ac.missions[imission].points[ip].engine

    if case == "design"
        #Design ducted fan for start of cruise
        ductedfancalc!(ac, case, imission, ip, initializes_engine)
        # ductedfancalc! dual-writes all outputs to both bare pare and typed state,
        # so eng_ip.design.{A2,A7,mbfD,pifD,NbfD} and eng_ip.{A7fac,Pfan} are current.

        parg[igA7] = eng_ip.design.A7 / eng_ip.A7fac

        pare[ieA2, :] .= pare[ieA2, ip]
        pare[ieA7, :] .= parg[igA7] .* pare[ieA7fac, :]

        pare[iembfD, :] .= pare[iembfD, ip]
        pare[iepifD, :] .= pare[iepifD, ip]
        pare[ieNbfD, :] .= pare[ieNbfD, ip]

        # Typed-state dual-writes: propagate design constants to all mission points.
        # Mirrors the bare pare broadcasts above (tasopt-048).
        for pt in ac.missions[imission].points
            pt.engine.design.A2   = eng_ip.design.A2
            pt.engine.design.A7   = parg[igA7] * pt.engine.A7fac
            pt.engine.design.mbfD = eng_ip.design.mbfD
            pt.engine.design.pifD = eng_ip.design.pifD
            pt.engine.design.NbfD = eng_ip.design.NbfD
        end

        #Design fuel cell for takeoff conditions
        ip_fcdes = iprotate
        eng_fcdes = ac.missions[imission].points[ip_fcdes].engine
        eng_fcdes.Pfanmax = pare[iePfanmax, ip_fcdes]   # sync Pfanmax to typed state
        Pfanmax = eng_fcdes.Pfanmax

        ## Model of electric machine to deliver Pfanmax
        Pmotormax = Pfanmax #100% efficiency for now
        ##
        fcdata.design_power = Pmotormax

        size_fuel_cell!(ac, ip_fcdes, imission)

        #Evaluate state at design point ip
        Pfan = eng_ip.Pfan
        ## Model of electric machine to deliver Pfan
        Pmotor = Pfan #100% efficiency for now
        ##
        fcdata.FC_power[ip, imission] = Pmotor
        operate_fuel_cell!(ac, ip, imission)

    else
        if ip in range(ipstatic, ipcutback)
            icall = 1

        elseif ip in range(ipclimb1, ipclimbn)
            icall = 2

            calculate_thrust_from_ROC!(ac, ip, imission) #Calculate thrust required for climb

        else
            icall = 2

        end

        ductedfancalc!(ac, case, imission, ip, initializes_engine)
        # ductedfancalc! dual-writes Pfan to both bare pare and typed state.

        Pfan = eng_ip.Pfan
        ## Model of electric machine to deliver Pfan
        Pmotor = Pfan #100% efficiency for now
        ##
        fcdata.FC_power[ip, imission] = Pmotor

        operate_fuel_cell!(ac, ip, imission)
    end
end

"""
    VerifyRadiatorHeat(engine, imission)

Verifies the heat balance between the fuel cell and the radiator for a given mission. This function 
compares the heat rejected by the fuel cell with the heat absorbed by the radiator at each mission 
point. If a significant imbalance is found, an error is raised.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `engine::Engine`: engine object containing fuel cell and heat exchanger data
    - `imission::Int64`: mission index

    **Output:**
    No direct outputs. If a heat imbalance is detected, the function throws an error with 
    information about the mission point where the imbalance occurred.

!!! note
    The comparison uses a relative tolerance (`rtol`) of `1e-6` to account for small numerical discrepancies.
"""
function VerifyRadiatorHeat(engine, imission)
    radiator = engine.heat_exchangers[1]
    fcdata = engine.data

    for ip in 1:iptotal
        if !isapprox(fcdata.FC_heat[ip, imission], radiator.HXgas_mission[ip, imission].Q, rtol=1e-6)
            error("Heat balance error in fuel cell & radiator")
        end
    end
end