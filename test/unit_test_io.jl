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
