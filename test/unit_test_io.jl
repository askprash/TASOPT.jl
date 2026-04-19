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
        @test saved["Propulsion"]["Turbomachinery"]["BPR"] ≈ ac_def.missions[1].points[1].engine.BPR
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
        eng_d = ac_def.missions[1].points[1].engine.design

        # Pressure ratios
        @test turb["diffuser_PR"]    ≈ eng_d.pid
        @test turb["burner_PR"]      ≈ eng_d.pib
        @test turb["fan_nozzle_PR"]  ≈ eng_d.pifn
        @test turb["core_nozzle_PR"] ≈ eng_d.pitn

        # Polytropic efficiencies
        @test turb["fan_eta_poly"]   ≈ eng_d.epolf
        @test turb["LPC_eta_poly"]   ≈ eng_d.epollc
        @test turb["HPC_eta_poly"]   ≈ eng_d.epolhc
        @test turb["HPT_eta_poly"]   ≈ eng_d.epolht
        @test turb["LPT_eta_poly"]   ≈ eng_d.epollt

        # Fan map constants
        @test turb["FPR0"]           ≈ eng_d.pifK
        @test turb["Kf_polyeff"]     ≈ eng_d.epfK

        # Duct Mach numbers
        @test turb["M2"]             ≈ eng_d.M2
        @test turb["M25"]            ≈ eng_d.M25

        # Spool losses
        @test turb["low_spool_loss"]  ≈ eng_d.epsl
        @test turb["high_spool_loss"] ≈ eng_d.epsh

        # Combustion efficiency
        @test comb["combustion_efficiency"] ≈ eng_d.etab

        # Cooling design parameters
        @test cool["hot_streak_T_allowance"] ≈ eng_d.dTstrk
        @test cool["M_turbine_blade_exit"]   ≈ eng_d.Mtexit
        @test cool["St"]                     ≈ eng_d.StA
        @test cool["e_film_cooling"]         ≈ eng_d.efilm
        @test cool["t_film_cooling"]         ≈ eng_d.tfilm

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

        let eng = ac_def.missions[1].points[1].engine
            @test turb["Fan_PR"]  ≈ eng.pif
            @test turb["LPC_PR"]  ≈ eng.pilc
            @test turb["OPR"]     ≈ eng.pilc * eng.pihc
            @test cool["M41"]               ≈ eng.design.M4a
            @test cool["cooling_air_V_ratio"] ≈ eng.design.ruc
            @test offt["Tt_offtake_air"]    ≈ eng.st9.Tt
            @test offt["Pt_offtake_air"]    ≈ eng.st9.pt
        end
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

        # Core nozzle area schedule — each flight point must match typed per-point A5fac
        @test cnoz["static"]       ≈ ac_def.missions[1].points[ipstatic].engine.A5fac
        @test cnoz["rotation"]     ≈ ac_def.missions[1].points[iprotate].engine.A5fac
        @test cnoz["cutback"]      ≈ ac_def.missions[1].points[ipcutback].engine.A5fac
        @test cnoz["climbstart"]   ≈ ac_def.missions[1].points[ipclimb1].engine.A5fac
        @test cnoz["climbend"]     ≈ ac_def.missions[1].points[ipclimbn].engine.A5fac
        @test cnoz["descentstart"] ≈ ac_def.missions[1].points[ipdescent1].engine.A5fac
        @test cnoz["descentend"]   ≈ ac_def.missions[1].points[ipdescentn].engine.A5fac

        # Fan nozzle area schedule — same 7 flight points
        @test fnoz["static"]       ≈ ac_def.missions[1].points[ipstatic].engine.A7fac
        @test fnoz["rotation"]     ≈ ac_def.missions[1].points[iprotate].engine.A7fac
        @test fnoz["cutback"]      ≈ ac_def.missions[1].points[ipcutback].engine.A7fac
        @test fnoz["climbstart"]   ≈ ac_def.missions[1].points[ipclimb1].engine.A7fac
        @test fnoz["climbend"]     ≈ ac_def.missions[1].points[ipclimbn].engine.A7fac
        @test fnoz["descentstart"] ≈ ac_def.missions[1].points[ipdescent1].engine.A7fac
        @test fnoz["descentend"]   ≈ ac_def.missions[1].points[ipdescentn].engine.A7fac

        # Tt4 at cruise and takeoff — saved as vectors over missions
        @test prop["Tt4_cruise"][1]  ≈ ac_def.missions[1].points[ipcruise1].engine.st4.Tt
        @test prop["Tt4_takeoff"][1] ≈ ac_def.missions[1].points[ipstatic].engine.st4.Tt

        rm(filepath_dw7)
    end

    # tasopt-1v7: verify that read_aircraft_model populates typed per-point engine
    # state for Tt4, T0 (takeoff), and nozzle area factors at input time.
    # After tasopt-j9l.45.14.4 bare pare is no longer populated at parse time;
    # only typed state carries these values until the first engine call.
    @testset "read_input populates per-point engine state from input (tasopt-1v7)" begin
        ac_1v7 = load_default_model()
        im = 1  # design mission

        # Tt4 at key flight points — typed-only after tasopt-j9l.45.14.4 (bare pare no longer
        # populated by read_input.jl; engine_state_to_pare! writes it on first engine call).
        # Tt4_takeoff = 1833.0 K; ipclimb1/ipcruise1 get Tt4_cruise = 1587.0 (read_input broadcast)
        @test ac_1v7.missions[im].points[ipstatic].engine.st4.Tt   ≈ 1833.0
        @test ac_1v7.missions[im].points[iptakeoff].engine.st4.Tt  ≈ 1833.0
        @test ac_1v7.missions[im].points[ipclimb1].engine.st4.Tt   ≈ 1587.0
        @test ac_1v7.missions[im].points[ipcruise1].engine.st4.Tt  ≈ 1587.0

        # T0 at takeoff points: 288.2 K (default_input.toml: takeoff_T = [288.2, 298.0])
        @test ac_1v7.missions[im].points[ipstatic].engine.T0  ≈ 288.2
        @test ac_1v7.missions[im].points[iprotate].engine.T0  ≈ 288.2
        @test ac_1v7.missions[im].points[iptakeoff].engine.T0 ≈ 288.2

        # Nozzle area factors: all 1.0 (default_input.toml uniform schedule)
        @test ac_1v7.missions[im].points[ipstatic].engine.A5fac   ≈ 1.0
        @test ac_1v7.missions[im].points[ipcruise1].engine.A5fac  ≈ 1.0
        @test ac_1v7.missions[im].points[ipclimb1].engine.A5fac   ≈ 1.0
        @test ac_1v7.missions[im].points[ipdescent1].engine.A5fac ≈ 1.0

        @test ac_1v7.missions[im].points[ipstatic].engine.A7fac   ≈ 1.0
        @test ac_1v7.missions[im].points[ipcruise1].engine.A7fac  ≈ 1.0
        @test ac_1v7.missions[im].points[ipclimb1].engine.A7fac   ≈ 1.0
        @test ac_1v7.missions[im].points[ipdescent1].engine.A7fac ≈ 1.0
    end

    # tasopt-50r: verify that read_aircraft_model populates typed design-point engine
    # state for turbomachinery scalars, cooling constants, and offtake discharge
    # conditions immediately after parsing — no sizing needed.
    @testset "read_input populates design-point engine state from input (tasopt-50r)" begin
        ac_50r = load_default_model()
        im = 1  # design mission
        eng_d = ac_50r.missions[im].points[ipcruise1].engine.design

        # Component pressure ratios — typed-only after tasopt-j9l.45.14.4
        # (bare pare no longer populated by read_input.jl)
        # Values from default_input.toml (same as tasopt-j9l.43)
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

        # Per-point pressure ratios and BPR — typed-only against TOML defaults
        # Fan_PR=1.685, LPC_PR=2.5, OPR=30 → pihc=30/2.5=12.0, BPR=5.1
        @test ac_50r.missions[im].points[ipcruise1].engine.pif  ≈ 1.685
        @test ac_50r.missions[im].points[ipcruise1].engine.pilc ≈ 2.5
        @test ac_50r.missions[im].points[ipcruise1].engine.pihc ≈ 12.0
        @test ac_50r.missions[im].points[ipcruise1].engine.BPR  ≈ 5.1
        @test ac_50r.missions[im].points[ipstatic].engine.BPR   ≈ 5.1

        # Cooling design constants
        @test eng_d.M4a    ≈ 0.9
        @test eng_d.ruc    ≈ 0.15
        @test eng_d.dTstrk ≈ 200.0   # hot_streak_T_allowance = 200.0 K
        @test eng_d.Mtexit ≈ 1.0
        @test eng_d.StA    ≈ 0.09
        @test eng_d.efilm  ≈ 0.70
        @test eng_d.tfilm  ≈ 0.30
        @test eng_d.fc0     ≈ 0.0
        @test eng_d.dehtdfc ≈ 0.0

        # Offtake discharge conditions (typed-only; values from default_input.toml)
        @test ac_50r.missions[im].points[ipcruise1].engine.st9.Tt ≈ 300.0
        @test ac_50r.missions[im].points[ipcruise1].engine.st9.pt ≈ 30e3
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

        # Typed-only after tasopt-j9l.45.14.4: bare pare no longer populated by read_input.jl.
        # All flight points carry fuel_temp = 280.0 K (default_input.toml)
        for ip in 1:iptotal
            @test ac_3ua.missions[im].points[ip].engine.Tfuel ≈ 280.0
        end

        # Round-trip check: save then reload preserves fuel_temp.
        filepath_3ua = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_3ua.toml")
        save_aircraft_model(ac_3ua, filepath_3ua)
        saved_toml = TOML.parsefile(filepath_3ua)
        @test saved_toml["Fuel"]["fuel_temp"] ≈ ac_3ua.missions[im].points[1].engine.Tfuel
        rm(filepath_3ua)
    end

    # tasopt-j9l.61/tasopt-w82: verify that read_input populates typed hvapcombustor at parse
    # time from the TOML fuel_enthalpy_vaporization value. Bare pare slot removed (tasopt-w82).
    @testset "read_input populates typed hvapcombustor from TOML (tasopt-j9l.61)" begin
        ac_j61 = load_default_model()
        im = 1  # design mission

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

    @test length(csv1[1]) == 784 # = entries w/ full ac `struct` and in default_output_indices
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
