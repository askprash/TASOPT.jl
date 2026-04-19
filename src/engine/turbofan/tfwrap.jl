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
    parg, _, para, pare, options, _, _, wing, _, _, engine = unpack_ac(ac, imission)
    
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

        pare_to_engine_state!(ac.missions[imission].points[ip].engine, view(pare, :, ip))
        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip), view(pare, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine_firstiter)

        # Propagate design-point scalars to all per-point typed EngineStates.
        # sync_design_scalars_to_pare! also writes the same 19 fields to bare
        # pare so that the pre-sync loop in mission_iteration.jl does not
        # overwrite mbfD/NbfD/etc. with zero values from an un-initialised bare
        # pare column (tasopt-j9l.45.14.2).  ruc/M4a are intentionally excluded
        # (written only for the design point inside tfcalc!; overwriting
        # non-design-point bare-pare columns triggers H2-sizing NaN).
        eng_ip = ac.missions[imission].points[ip].engine
        parg[igA5] = eng_ip.design.A5 / eng_ip.A5fac
        parg[igA7] = eng_ip.design.A7 / eng_ip.A7fac

        for jp = 1:iptotal
            eng_jp = ac.missions[imission].points[jp].engine
            eng_jp.design.A2    = eng_ip.design.A2
            eng_jp.design.A25   = eng_ip.design.A25
            eng_jp.design.A5    = parg[igA5] * eng_jp.A5fac
            eng_jp.design.A7    = parg[igA7] * eng_jp.A7fac
            eng_jp.design.NbfD  = eng_ip.design.NbfD
            eng_jp.design.NblcD = eng_ip.design.NblcD
            eng_jp.design.NbhcD = eng_ip.design.NbhcD
            eng_jp.design.NbhtD = eng_ip.design.NbhtD
            eng_jp.design.NbltD = eng_ip.design.NbltD
            eng_jp.design.mbfD  = eng_ip.design.mbfD
            eng_jp.design.mblcD = eng_ip.design.mblcD
            eng_jp.design.mbhcD = eng_ip.design.mbhcD
            eng_jp.design.mbhtD = eng_ip.design.mbhtD
            eng_jp.design.mbltD = eng_ip.design.mbltD
            eng_jp.design.pifD  = eng_ip.design.pifD
            eng_jp.design.pilcD = eng_ip.design.pilcD
            eng_jp.design.pihcD = eng_ip.design.pihcD
            eng_jp.design.pihtD = eng_ip.design.pihtD
            eng_jp.design.piltD = eng_ip.design.piltD
            sync_design_scalars_to_pare!(eng_jp.design, view(pare, :, jp))
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

        pare_to_engine_state!(ac.missions[imission].points[ip].engine, view(pare, :, ip))
        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip), view(pare, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine)

    elseif case == "cooling_sizing"
        opt_calc_call = CalcMode.FixedTt4OffDes
        opt_cooling = CoolingOpt.FixedTmetal
        pare_to_engine_state!(ac.missions[imission].points[ip].engine, view(pare, :, ip))
        ichoke5, ichoke7 = tfcalc!(wing, engine, parg, view(para, :, ip), view(pare, :, ip),
            ac.missions[imission].points[ip].engine,
            ip, options.ifuel, opt_calc_call, opt_cooling, initializes_engine)

        # Tmetal was specified... propagate blade-row cooling fractions to all points.
        # sync_cooling_scalars_to_pare! keeps bare pare in sync so that the
        # pre-sync in mission_iteration.jl does not clobber the freshly computed
        # values with stale bare-pare values (tasopt-j9l.45.14.2).
        eng_ip = ac.missions[imission].points[ip].engine
        for jp = 1:iptotal
            eng_jp = ac.missions[imission].points[jp].engine
            eng_jp.design.epsrow = eng_ip.design.epsrow
            eng_jp.design.fc     = eng_ip.design.fc
            sync_cooling_scalars_to_pare!(eng_jp.design, view(pare, :, jp))
        end
    end
end