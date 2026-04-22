""" 
    tfwrap!(ac, case, imission, ip, initializes_engine, iterw = 0)

Calls the turbofan sizing or off-design performance functions for the aircraft's
turbofan model. This function is basically a wrapper on tfcalc!, going from the
basic engine inputs to those required by the function and storing the outputs.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `case::String`: case identifier, e.g. "sizing" or "off_design"
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `initializes_engine::Bool`: flag to initialize engine:
       - `true`: initialize variables for iteration in engine
       - `false`: use current variables as initial guesses in engine
    - `iterw::Int64`: sizing loop iteration

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function tfwrap!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)
    #Unpack data storage arrays
    parg, _, para, options, _, _, wing, _, _, engine, _ = unpack_ac(ac, imission)
    
    if case == "design"
        opt_calc_call = CalcMode.Sizing
        opt_cooling = CoolingOpt.FixedCoolingFlowRatio
        if (iterw == 1 || (initializes_engine))
            # initialize engine state
            initializes_engine_firstiter  = true
        else
            # start with current engine state
            initializes_engine_firstiter  = false
        end

        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine_firstiter)

        # Propagate design-point scalars to all per-point typed EngineStates.
        # ruc/M4a are intentionally excluded (written only for the design point
        # inside tfcalc!).
        eng_ip = ac.missions[imission].points[ip].engine
        parg[igA5] = eng_ip.design.A8 / eng_ip.A5fac
        parg[igA7] = eng_ip.design.A18 / eng_ip.A7fac

        for jp = 1:iptotal
            eng_jp = ac.missions[imission].points[jp].engine
            eng_jp.design.A2    = eng_ip.design.A2
            eng_jp.design.A25   = eng_ip.design.A25
            eng_jp.design.A8    = parg[igA5] * eng_jp.A5fac
            eng_jp.design.A18   = parg[igA7] * eng_jp.A7fac
            eng_jp.design.Nb_fan_des  = eng_ip.design.Nb_fan_des
            eng_jp.design.Nb_lpc_des = eng_ip.design.Nb_lpc_des
            eng_jp.design.Nb_hpc_des = eng_ip.design.Nb_hpc_des
            eng_jp.design.Nb_hpt_des = eng_ip.design.Nb_hpt_des
            eng_jp.design.Nb_lpt_des = eng_ip.design.Nb_lpt_des
            eng_jp.design.mb_fan_des  = eng_ip.design.mb_fan_des
            eng_jp.design.mb_lpc_des = eng_ip.design.mb_lpc_des
            eng_jp.design.mb_hpc_des = eng_ip.design.mb_hpc_des
            eng_jp.design.mb_hpt_des = eng_ip.design.mb_hpt_des
            eng_jp.design.mb_lpt_des = eng_ip.design.mb_lpt_des
            eng_jp.design.pi_fan_des  = eng_ip.design.pi_fan_des
            eng_jp.design.pi_lpc_des = eng_ip.design.pi_lpc_des
            eng_jp.design.pi_hpc_des = eng_ip.design.pi_hpc_des
            eng_jp.design.pi_hpt_des = eng_ip.design.pi_hpt_des
            eng_jp.design.pi_lpt_des = eng_ip.design.pi_lpt_des
        end
        
    elseif case == "off_design"
        #assume operating at max allowable temp if during TO and climb
        if ip in range(ipstatic, ipclimbn)
            opt_calc_call = CalcMode.FixedTt4OffDes
        #otherwise, thrust balance sets op point
        else
            opt_calc_call = CalcMode.FixedFeOffDes
        end
        opt_cooling = CoolingOpt.FixedCoolingFlowRatio

        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine)

    elseif case == "cooling_sizing"
        opt_calc_call = CalcMode.FixedTt4OffDes
        opt_cooling = CoolingOpt.FixedTmetal
        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine)

        # Tmetal was specified... propagate blade-row cooling fractions to all points.
        eng_ip = ac.missions[imission].points[ip].engine
        for jp = 1:iptotal
            eng_jp = ac.missions[imission].points[jp].engine
            eng_jp.design.epsrow = eng_ip.design.epsrow
            eng_jp.design.fc     = eng_ip.design.fc
        end
    end
end