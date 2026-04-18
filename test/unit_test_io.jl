#tests io functionalities
@testset "io" verbose=true begin
#A: readable TOML saves
    #check that the default model is sized identically via MTOW
        #when round-tripped via model save and read
    ac_def = load_default_model()

    # tasopt-rcy: verify BPR is read from typed engine state in save_aircraft_model.
    # The sync in save_aircraft_model populates typed state from pare before reading,
    # so the saved BPR must equal the pare design value whether or not sized.
    @testset "save_model reads BPR from typed engine state (tasopt-rcy)" begin
        import TOML
        filepath_bpr = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_bpr.toml")
        save_aircraft_model(ac_def, filepath_bpr)
        saved = TOML.parsefile(filepath_bpr)
        @test saved["Propulsion"]["Turbomachinery"]["BPR"] ≈ ac_def.pare[ieBPR, 1, 1]
        rm(filepath_bpr)
    end

    # tasopt-86a: verify the 21 design constants newly migrated from bare pare
    # reads to eng.design.* are saved correctly.
    @testset "save_model reads expanded design constants from typed state (tasopt-86a)" begin
        import TOML
        filepath_86a = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_86a.toml")
        save_aircraft_model(ac_def, filepath_86a)
        saved = TOML.parsefile(filepath_86a)
        turb = saved["Propulsion"]["Turbomachinery"]
        comb = saved["Propulsion"]["Combustor"]
        cool = saved["Propulsion"]["Cooling"]

        # Pressure ratios
        @test turb["diffuser_PR"]    ≈ ac_def.pare[iepid,    1, 1]
        @test turb["burner_PR"]      ≈ ac_def.pare[iepib,    1, 1]
        @test turb["fan_nozzle_PR"]  ≈ ac_def.pare[iepifn,   1, 1]
        @test turb["core_nozzle_PR"] ≈ ac_def.pare[iepitn,   1, 1]

        # Polytropic efficiencies
        @test turb["fan_eta_poly"]   ≈ ac_def.pare[ieepolf,  1, 1]
        @test turb["LPC_eta_poly"]   ≈ ac_def.pare[ieepollc, 1, 1]
        @test turb["HPC_eta_poly"]   ≈ ac_def.pare[ieepolhc, 1, 1]
        @test turb["HPT_eta_poly"]   ≈ ac_def.pare[ieepolht, 1, 1]
        @test turb["LPT_eta_poly"]   ≈ ac_def.pare[ieepollt, 1, 1]

        # Fan map constants
        @test turb["FPR0"]           ≈ ac_def.pare[iepifK,   1, 1]
        @test turb["Kf_polyeff"]     ≈ ac_def.pare[ieepfK,   1, 1]

        # Duct Mach numbers
        @test turb["M2"]             ≈ ac_def.pare[ieM2,     1, 1]
        @test turb["M25"]            ≈ ac_def.pare[ieM25,    1, 1]

        # Spool losses
        @test turb["low_spool_loss"]  ≈ ac_def.pare[ieepsl,  1, 1]
        @test turb["high_spool_loss"] ≈ ac_def.pare[ieepsh,  1, 1]

        # Combustion efficiency
        @test comb["combustion_efficiency"] ≈ ac_def.pare[ieetab, 1, 1]

        # Cooling design parameters
        @test cool["hot_streak_T_allowance"] ≈ ac_def.pare[iedTstrk, 1, 1]
        @test cool["M_turbine_blade_exit"]   ≈ ac_def.pare[ieMtexit, 1, 1]
        @test cool["St"]                     ≈ ac_def.pare[ieStA,    1, 1]
        @test cool["e_film_cooling"]         ≈ ac_def.pare[ieefilm,  1, 1]
        @test cool["t_film_cooling"]         ≈ ac_def.pare[ietfilm,  1, 1]

        rm(filepath_86a)
    end

    # tasopt-sxv: verify design-point scalars (Fan_PR, LPC_PR, OPR, M41,
    # cooling_air_V_ratio, Tt/Pt offtake) are read from typed engine state.
    @testset "save_model reads design scalars from typed engine state (tasopt-sxv)" begin
        import TOML
        filepath_sxv = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_sxv.toml")
        save_aircraft_model(ac_def, filepath_sxv)
        saved = TOML.parsefile(filepath_sxv)
        turb = saved["Propulsion"]["Turbomachinery"]
        cool = saved["Propulsion"]["Cooling"]
        offt = saved["Propulsion"]["Offtakes"]

        @test turb["Fan_PR"]  ≈ ac_def.pare[iepif,  1, 1]
        @test turb["LPC_PR"]  ≈ ac_def.pare[iepilc, 1, 1]
        @test turb["OPR"]     ≈ ac_def.pare[iepilc, 1, 1] * ac_def.pare[iepihc, 1, 1]
        @test cool["M41"]               ≈ ac_def.pare[ieM4a,  1, 1]
        @test cool["cooling_air_V_ratio"] ≈ ac_def.pare[ieruc,   1, 1]
        @test offt["Tt_offtake_air"]    ≈ ac_def.pare[ieTt9,  1, 1]
        @test offt["Pt_offtake_air"]    ≈ ac_def.pare[iept9,  1, 1]
        rm(filepath_sxv)
    end

    # tasopt-dw7: verify per-point nozzle area factors and Tt4 are saved from
    # typed per-point engine state rather than bare pare reads.
    @testset "save_model reads nozzle schedule and Tt4 from typed state (tasopt-dw7)" begin
        import TOML
        filepath_dw7 = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_dw7.toml")
        save_aircraft_model(ac_def, filepath_dw7)
        saved = TOML.parsefile(filepath_dw7)
        prop  = saved["Propulsion"]
        cnoz  = prop["Nozzles"]["core_nozzle_area"]
        fnoz  = prop["Nozzles"]["fan_nozzle_area"]

        # Core nozzle area schedule — each flight point must match pare[ieA5fac, ip, 1]
        @test cnoz["static"]       ≈ ac_def.pare[ieA5fac, ipstatic,    1]
        @test cnoz["rotation"]     ≈ ac_def.pare[ieA5fac, iprotate,    1]
        @test cnoz["cutback"]      ≈ ac_def.pare[ieA5fac, ipcutback,   1]
        @test cnoz["climbstart"]   ≈ ac_def.pare[ieA5fac, ipclimb1,    1]
        @test cnoz["climbend"]     ≈ ac_def.pare[ieA5fac, ipclimbn,    1]
        @test cnoz["descentstart"] ≈ ac_def.pare[ieA5fac, ipdescent1,  1]
        @test cnoz["descentend"]   ≈ ac_def.pare[ieA5fac, ipdescentn,  1]

        # Fan nozzle area schedule — same 7 flight points
        @test fnoz["static"]       ≈ ac_def.pare[ieA7fac, ipstatic,    1]
        @test fnoz["rotation"]     ≈ ac_def.pare[ieA7fac, iprotate,    1]
        @test fnoz["cutback"]      ≈ ac_def.pare[ieA7fac, ipcutback,   1]
        @test fnoz["climbstart"]   ≈ ac_def.pare[ieA7fac, ipclimb1,    1]
        @test fnoz["climbend"]     ≈ ac_def.pare[ieA7fac, ipclimbn,    1]
        @test fnoz["descentstart"] ≈ ac_def.pare[ieA7fac, ipdescent1,  1]
        @test fnoz["descentend"]   ≈ ac_def.pare[ieA7fac, ipdescentn,  1]

        # Tt4 at cruise and takeoff — saved as vectors over missions
        @test prop["Tt4_cruise"][1]  ≈ ac_def.pare[ieTt4, ipcruise1, 1]
        @test prop["Tt4_takeoff"][1] ≈ ac_def.pare[ieTt4, ipstatic,  1]

        rm(filepath_dw7)
    end

    # tasopt-1v7: verify that read_aircraft_model populates typed per-point engine
    # state for Tt4, T0 (takeoff), and nozzle area factors at input time.
    # Typed state must agree with pare immediately after parsing — no sizing needed.
    @testset "read_input populates per-point engine state from input (tasopt-1v7)" begin
        ac_1v7 = load_default_model()
        im = 1  # design mission

        # Tt4 at all key flight points: cruise, takeoff
        @test ac_1v7.missions[im].points[ipcruise1].engine.st4.Tt ≈ ac_1v7.pare[ieTt4, ipcruise1, im]
        @test ac_1v7.missions[im].points[ipstatic].engine.st4.Tt  ≈ ac_1v7.pare[ieTt4, ipstatic,  im]
        @test ac_1v7.missions[im].points[iptakeoff].engine.st4.Tt ≈ ac_1v7.pare[ieTt4, iptakeoff, im]
        @test ac_1v7.missions[im].points[ipclimb1].engine.st4.Tt  ≈ ac_1v7.pare[ieTt4, ipclimb1,  im]

        # T0 at takeoff flight points
        @test ac_1v7.missions[im].points[ipstatic].engine.T0  ≈ ac_1v7.pare[ieT0, ipstatic,  im]
        @test ac_1v7.missions[im].points[iprotate].engine.T0  ≈ ac_1v7.pare[ieT0, iprotate,  im]
        @test ac_1v7.missions[im].points[iptakeoff].engine.T0 ≈ ac_1v7.pare[ieT0, iptakeoff, im]

        # Nozzle area factors at representative flight points
        @test ac_1v7.missions[im].points[ipstatic].engine.A5fac   ≈ ac_1v7.pare[ieA5fac, ipstatic,   im]
        @test ac_1v7.missions[im].points[ipcruise1].engine.A5fac  ≈ ac_1v7.pare[ieA5fac, ipcruise1,  im]
        @test ac_1v7.missions[im].points[ipclimb1].engine.A5fac   ≈ ac_1v7.pare[ieA5fac, ipclimb1,   im]
        @test ac_1v7.missions[im].points[ipdescent1].engine.A5fac ≈ ac_1v7.pare[ieA5fac, ipdescent1, im]

        @test ac_1v7.missions[im].points[ipstatic].engine.A7fac   ≈ ac_1v7.pare[ieA7fac, ipstatic,   im]
        @test ac_1v7.missions[im].points[ipcruise1].engine.A7fac  ≈ ac_1v7.pare[ieA7fac, ipcruise1,  im]
        @test ac_1v7.missions[im].points[ipclimb1].engine.A7fac   ≈ ac_1v7.pare[ieA7fac, ipclimb1,   im]
        @test ac_1v7.missions[im].points[ipdescent1].engine.A7fac ≈ ac_1v7.pare[ieA7fac, ipdescent1, im]
    end

    # tasopt-50r: verify that read_aircraft_model populates typed design-point engine
    # state for turbomachinery scalars, cooling constants, and offtake discharge
    # conditions immediately after parsing — no sizing needed.
    @testset "read_input populates design-point engine state from input (tasopt-50r)" begin
        ac_50r = load_default_model()
        im = 1  # design mission
        eng_d = ac_50r.missions[im].points[ipcruise1].engine.design

        # Component pressure ratios
        @test eng_d.pid    ≈ ac_50r.pare[iepid,    1, im]
        @test eng_d.pib    ≈ ac_50r.pare[iepib,    1, im]
        @test eng_d.pifn   ≈ ac_50r.pare[iepifn,   1, im]
        @test eng_d.pitn   ≈ ac_50r.pare[iepitn,   1, im]

        # Polytropic efficiencies
        @test eng_d.epolf  ≈ ac_50r.pare[ieepolf,  1, im]
        @test eng_d.epollc ≈ ac_50r.pare[ieepollc, 1, im]
        @test eng_d.epolhc ≈ ac_50r.pare[ieepolhc, 1, im]
        @test eng_d.epolht ≈ ac_50r.pare[ieepolht, 1, im]
        @test eng_d.epollt ≈ ac_50r.pare[ieepollt, 1, im]

        # Combustion efficiency, duct Mach numbers, spool losses
        @test eng_d.etab   ≈ ac_50r.pare[ieetab,   1, im]
        @test eng_d.M2     ≈ ac_50r.pare[ieM2,     1, im]
        @test eng_d.M25    ≈ ac_50r.pare[ieM25,    1, im]
        @test eng_d.epsl   ≈ ac_50r.pare[ieepsl,   1, im]
        @test eng_d.epsh   ≈ ac_50r.pare[ieepsh,   1, im]

        # Per-point pressure ratios and BPR (uniform at parse time)
        @test ac_50r.missions[im].points[ipcruise1].engine.pif  ≈ ac_50r.pare[iepif,  ipcruise1, im]
        @test ac_50r.missions[im].points[ipcruise1].engine.pilc ≈ ac_50r.pare[iepilc, ipcruise1, im]
        @test ac_50r.missions[im].points[ipcruise1].engine.pihc ≈ ac_50r.pare[iepihc, ipcruise1, im]
        @test ac_50r.missions[im].points[ipcruise1].engine.BPR  ≈ ac_50r.pare[ieBPR,  ipcruise1, im]
        @test ac_50r.missions[im].points[ipstatic].engine.BPR   ≈ ac_50r.pare[ieBPR,  ipstatic,  im]

        # Cooling design constants
        @test eng_d.M4a    ≈ ac_50r.pare[ieM4a,    1, im]
        @test eng_d.ruc    ≈ ac_50r.pare[ieruc,    1, im]
        @test eng_d.dTstrk ≈ ac_50r.pare[iedTstrk, 1, im]
        @test eng_d.Mtexit ≈ ac_50r.pare[ieMtexit, 1, im]
        @test eng_d.StA    ≈ ac_50r.pare[ieStA,    1, im]
        @test eng_d.efilm  ≈ ac_50r.pare[ieefilm,  1, im]
        @test eng_d.tfilm  ≈ ac_50r.pare[ietfilm,  1, im]
        @test eng_d.fc0     ≈ ac_50r.pare[iefc0,     1, im]
        @test eng_d.dehtdfc ≈ ac_50r.pare[iedehtdfc, 1, im]

        # Offtake discharge conditions
        @test ac_50r.missions[im].points[ipcruise1].engine.st9.Tt ≈ ac_50r.pare[ieTt9, 1, im]
        @test ac_50r.missions[im].points[ipcruise1].engine.st9.pt ≈ ac_50r.pare[iept9, 1, im]
    end

    # tasopt-j9l.43: verify that read_input directly populates typed design-point engine
    # state, and that the values match the known defaults from default_input.toml.
    # This pins the data-flow: typed state must be populated at parse time from the TOML
    # values, not merely be consistent with pare (which the tasopt-50r test already checks).
    @testset "read_input direct-populates design-point engine state (tasopt-j9l.43)" begin
        ac_j43 = load_default_model()
        im = 1
        eng_d = ac_j43.missions[im].points[ipcruise1].engine.design

        # Turbomachinery component pressure ratios (default_input.toml values)
        @test eng_d.pid    ≈ 0.998
        @test eng_d.pib    ≈ 0.94
        @test eng_d.pifn   ≈ 0.98
        @test eng_d.pitn   ≈ 0.989

        # Polytropic efficiencies
        @test eng_d.epolf  ≈ 0.8948
        @test eng_d.epollc ≈ 0.88
        @test eng_d.epolhc ≈ 0.87
        @test eng_d.epolht ≈ 0.889
        @test eng_d.epollt ≈ 0.899

        # Combustion efficiency, duct Mach numbers, spool losses
        @test eng_d.etab   ≈ 0.98
        @test eng_d.M2     ≈ 0.60
        @test eng_d.M25    ≈ 0.60
        @test eng_d.epsl   ≈ 0.01
        @test eng_d.epsh   ≈ 0.022

        # Cooling design constants
        @test eng_d.M4a    ≈ 0.9
        @test eng_d.ruc    ≈ 0.15
        @test eng_d.Mtexit ≈ 1.0
        @test eng_d.StA    ≈ 0.09
        @test eng_d.efilm  ≈ 0.70
        @test eng_d.tfilm  ≈ 0.30
        @test eng_d.fc0     ≈ 0.0
        @test eng_d.dehtdfc ≈ 0.0

        # Offtake discharge conditions (typed st9)
        eng_pt = ac_j43.missions[im].points[ipcruise1].engine
        @test eng_pt.st9.Tt ≈ 300.0
        @test eng_pt.st9.pt ≈ 30e3
    end

    # tasopt-3ua: verify that read_aircraft_model mirrors fuel temperature to typed
    # engine state for all flight points, and that save_model reads it from typed state.
    @testset "read_input and save_model use typed Tfuel (tasopt-3ua)" begin
        ac_3ua = load_default_model()
        im = 1  # design mission

        # Mirror check: typed Tfuel equals bare pare at every flight point after parse.
        for ip in 1:iptotal
            @test ac_3ua.missions[im].points[ip].engine.Tfuel ≈ ac_3ua.pare[ieTfuel, ip, im]
        end

        # Round-trip check: save then reload preserves fuel_temp.
        filepath_3ua = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_3ua.toml")
        save_aircraft_model(ac_3ua, filepath_3ua)
        saved_toml = TOML.parsefile(filepath_3ua)
        @test saved_toml["Fuel"]["fuel_temp"] ≈ ac_3ua.pare[ieTfuel, 1, im]
        rm(filepath_3ua)
    end

    # tasopt-j9l.61: verify that read_input populates typed hvapcombustor at parse time
    # from the TOML fuel_enthalpy_vaporization value, so that typed state is consistent
    # with pare[iehvapcombustor] immediately after parsing (before any HX or sizing call).
    @testset "read_input populates typed hvapcombustor from TOML (tasopt-j9l.61)" begin
        ac_j61 = load_default_model()
        im = 1  # design mission

        # Agreement with pare: typed hvapcombustor must equal bare pare at every flight point.
        for ip in 1:iptotal
            @test ac_j61.missions[im].points[ip].engine.hvapcombustor ≈
                  ac_j61.pare[iehvapcombustor, ip, im]
        end

        # Absolute value: default_input.toml has fuel_enthalpy_vaporization = 0.0 J/kg.
        # This pins that the default flows through correctly.
        @test ac_j61.missions[im].points[ipcruise1].engine.hvapcombustor ≈ 0.0
    end

    # tasopt-j9l.44: verify that save_model writes the constant_TSFC schedule (climb/cruise/descent
    # TSFC and rate_of_climb) from typed engine state when the engine architecture is ConstantTSFC.
    # The round-trip must preserve all four TOML keys under [Propulsion].
    @testset "save_model writes constant_TSFC schedule from typed state (tasopt-j9l.44)" begin
        import TOML
        ac_ctsfc = load_default_model()
        im = 1  # design mission

        # Patch architecture to ConstantTSFC.
        ac_ctsfc.options.opt_prop_sys_arch = PropSysArch.ConstantTSFC

        # Inject known TSFC schedule into typed state AND pare per-point.
        # save_aircraft_model calls pare_to_engine_state! for several points
        # (nozzle sync, Tt4 sync) which overwrites typed state from pare; so both
        # must agree for the assertion to hold after save.
        tsfc_climb   = 1.8e-4  # [kg/N/s = 1/s for consistent units]
        tsfc_cruise  = 1.6e-4
        tsfc_descent = 2.1e-4
        roc_ms       = 500 * 0.3048 / 60  # 500 ft/min → m/s

        for ip in ipclimb1:ipclimbn
            ac_ctsfc.missions[im].points[ip].engine.TSFC = tsfc_climb
            ac_ctsfc.pare[ieTSFC, ip, im] = tsfc_climb
        end
        for ip in ipcruise1:ipcruisen
            ac_ctsfc.missions[im].points[ip].engine.TSFC = tsfc_cruise
            ac_ctsfc.pare[ieTSFC, ip, im] = tsfc_cruise
        end
        for ip in ipdescent1:ipdescentn
            ac_ctsfc.missions[im].points[ip].engine.TSFC = tsfc_descent
            ac_ctsfc.pare[ieTSFC, ip, im] = tsfc_descent
        end
        ac_ctsfc.para[iaROCdes, ipclimb1:ipclimbn, im] .= roc_ms

        # Save and parse the TOML.
        filepath_j44 = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_j9l44.toml")
        save_aircraft_model(ac_ctsfc, filepath_j44)
        saved = TOML.parsefile(filepath_j44)
        prop  = saved["Propulsion"]

        @test prop["climb_TSFC"]   ≈ tsfc_climb
        @test prop["cruise_TSFC"]  ≈ tsfc_cruise
        @test prop["descent_TSFC"] ≈ tsfc_descent
        @test prop["rate_of_climb"] ≈ roc_ms

        # Turbofan default must NOT write TSFC schedule keys.
        filepath_j44_tf = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_j9l44_tf.toml")
        save_aircraft_model(load_default_model(), filepath_j44_tf)
        saved_tf = TOML.parsefile(filepath_j44_tf)
        @test !haskey(saved_tf["Propulsion"], "climb_TSFC")
        @test !haskey(saved_tf["Propulsion"], "cruise_TSFC")
        @test !haskey(saved_tf["Propulsion"], "descent_TSFC")

        rm(filepath_j44)
        rm(filepath_j44_tf)
    end

    # tasopt-j9l.62: verify that read_input populates Gearf/HTRf/HTRlc/HTRhc into typed
    # engine DesignState, and that save_aircraft_model reads them from there (not parg).
    @testset "turbomachinery geometry in typed DesignState (tasopt-j9l.62)" begin
        import TOML

        # --- read side: typed fields match known default_input.toml constants ---
        eng_d = ac_def.missions[1].points[1].engine
        @test eng_d.design.Gearf  ≈ 1.0   # default: direct-drive
        @test eng_d.design.HTRf   ≈ 0.30
        @test eng_d.design.HTRlc  ≈ 0.60
        @test eng_d.design.HTRhc  ≈ 0.80

        # uniform across all flight points and missions
        for im in 1:size(ac_def.pare, 3)
            for ip in 1:size(ac_def.pare, 2)
                e = ac_def.missions[im].points[ip].engine
                @test e.design.Gearf  == eng_d.design.Gearf
                @test e.design.HTRf   == eng_d.design.HTRf
                @test e.design.HTRlc  == eng_d.design.HTRlc
                @test e.design.HTRhc  == eng_d.design.HTRhc
            end
        end

        # --- write side: save_aircraft_model writes from typed state ---
        filepath_j62 = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_j9l62.toml")
        save_aircraft_model(ac_def, filepath_j62)
        saved_j62 = TOML.parsefile(filepath_j62)
        pt = saved_j62["Propulsion"]["Turbomachinery"]
        @test pt["gear_ratio"] ≈ eng_d.design.Gearf
        @test pt["HTR_fan"]   ≈ eng_d.design.HTRf
        @test pt["HTR_LPC"]   ≈ eng_d.design.HTRlc
        @test pt["HTR_HPC"]   ≈ eng_d.design.HTRhc
        rm(filepath_j62)
    end

    filepath_rewrite = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_rewrite.toml")
    save_aircraft_model(ac_def, filepath_rewrite)
    ac_reread = read_aircraft_model(filepath_rewrite)

    size_aircraft!(ac_def, Ldebug=false, printiter=false)
    size_aircraft!(ac_reread, Ldebug=false, printiter=false)

    @test ac_def.parg[igWMTO] ≈ ac_reread.parg[igWMTO]
    rm(filepath_rewrite)

    #check via MTOW that changing an important parameter survives the save
    # and changes the solution
    ac_lopay = load_default_model()
    ac_lopay.parm[imWpay] = ac_def.parm[imWperpax, 1] * 30 #thirty pax in N
    filepath_nopay = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_nopay.toml")
    save_aircraft_model(ac_lopay, filepath_nopay)

    ac_lopay_reread = read_aircraft_model(filepath_nopay)
    size_aircraft!(ac_lopay, Ldebug=false, printiter=false)
    size_aircraft!(ac_lopay_reread, Ldebug=false, printiter=false)
    
    @test ac_lopay.parg[igWMTO] ≈ ac_lopay_reread.parg[igWMTO]
    @test !(ac_lopay.parg[igWMTO] ≈ ac_def.parg[igWMTO])
    rm(filepath_nopay)

#B: quicksaves and loads
    #test that quicksave/load roundtrip default aircraft sizes identically to default load
    filepath_quick = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_quick.toml")
    quicksave_aircraft(load_default_model(), filepath_quick)
    ac_quick = quickload_aircraft(filepath_quick)
    size_aircraft!(ac_quick, Ldebug=false, printiter=false)

    @test ac_quick.parg[igWMTO] ≈ ac_def.parg[igWMTO]
    rm(filepath_quick)

    #check via MTOW that changing an important parameter survives the quicksave
    # and changes the solution
    filepath_quick_nopay = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_quick_nopay.toml")
    ac_quick.parm[imWpay] = ac_def.parm[imWperpax, 1] * 30 #thirty pax in N
    quicksave_aircraft(ac_quick, filepath_quick_nopay)

    ac_quick_nopay_reread = quickload_aircraft(filepath_quick_nopay)
    size_aircraft!(ac_quick_nopay_reread, Ldebug=false, printiter=false)
    @test ac_quick_nopay_reread.parg[igWMTO] ≈ ac_lopay.parg[igWMTO]
    rm(filepath_quick_nopay)

#C: outputs to .csv
    using CSV
    #generate file paths
    filepath_csv = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_def.csv")
    suffix = "_1"
    insertpos = length(filepath_csv) - 3
    filepath_csv2 = filepath_csv[1:insertpos-1] * suffix * filepath_csv[insertpos:end]

    #cleanup in case test was interrupted before
    isfile(filepath_csv) ? rm(filepath_csv) : nothing
    isfile(filepath_csv2) ? rm(filepath_csv2) : nothing
    
    #generate output files, indices are default_output_indices
    output_csv(ac_def, filepath_csv)
    output_csv(ac_def, filepath_csv, 
                includeMissions = true, includeFlightPoints = true)
    output_csv(ac_def, filepath_csv,
                includeMissions = false, includeFlightPoints = true,
                forceMatrices = true)
    output_csv(ac_def, filepath_csv, includeFlightPoints = true)

    #this call generates the second file since indices don't match and fuse_tank excluded
    output_csv(ac_def, filepath_csv, indices = Dict(), struct_excludes = ["fuse_tank"])
    
    #test that it creates the files as expected
    @test isfile(filepath_csv)
    @test isfile(filepath_csv2)
    
    #pull the files
    csv1 = CSV.File(filepath_csv)
    csv2 = CSV.File(filepath_csv2)

    #check row and column counts
    @test size(csv1,1) == 4 #4 rows w default indices
    @test size(csv2,1) == 1 #1 row with addl indices (all)

    @test length(csv1[1]) == 785 # = entries w/ full ac `struct` and in default_output_indices
    @test length(csv2[1]) == 1155 # = entries w/ full ac `struct` and all output_indices

    #test the nested vectors within par arrays
    #a: row 1 in both csvs matches the design cruise point/mission 
    @test parse(Float64, string(csv1[1].iaalt)) == ac_def.para[iaalt,ipcruise1,1]
    @test parse(Float64, string(csv2[1].iaalt)) == ac_def.para[iaalt,ipcruise1,1]
    
    #b: row 2 has the correct structure (m flight points, n missions)
    #note - for simplicity of imports, evaluate structure by counting brackets
    # since much of the data is parsed as Strings when using CSV.File
    #if more than one mission, one more bracket added
    entry2 = csv1[2].iaalt
    bracket_for_missions = Int( (size(ac_def.parm,2) > 1) ) 
    @test count(ichar->(ichar=='['), entry2) == iptotal + bracket_for_missions
    #c: row 3 has the same structure as b despite 1 mission bc forceMatrices
    entry3 = csv1[3].iaalt
    @test count(ichar->(ichar=='['), entry3) == iptotal + bracket_for_missions
    #d: row 4 only has flight points, thus 1 bracket
    entry4 = csv1[4].iaalt
    @test count(ichar->(ichar=='['), entry4) == 1

    #test that structs are properly saved
    #a: check fuselage detail, wing detail
    @test csv1[1][Symbol("fuselage.layout.cross_section.radius")] == ac_def.fuselage.layout.cross_section.radius
    @test csv2[1][Symbol("wing.inboard.cross_section.thickness_to_chord")] == ac_def.wing.inboard.cross_section.thickness_to_chord

    #b: test exclusion of structs when indicated
    @test !haskey(csv2[1], Symbol("fuse_tank.tank_type")) #excluded

    #cleanup
    rm(filepath_csv)
    rm(filepath_csv2)

end #testset "io"
