# Define a function to check if each value in two structs is equal
function check_struct_equivalence(s1, s2)
    fields_s1 = fieldnames(typeof(s1))
    fields_s2 = fieldnames(typeof(s2))
    
    # Check if both structs have the same fields
    if fields_s1 != fields_s2
        return false
    end
    
    # Check if each field has the same value in both structs
    for field in fields_s1
        val1 = getproperty(s1, field)
        val2 = getproperty(s2, field)
        if typeof(val1) == typeof(val2)
            if typeof(val1) != Float64
                if !check_struct_equivalence(val1, val2)
                    return false
                end
            else
                # println(field)
                @test val1 ≈ val2 
            end
        else
            return false
        end
    end
    
    return true
end

#Simple function to call fly_mission!() and test on- and off-design performance
function test_ac_off_design(ac, PFEI, Wfuel, WTO)
    @testset "Off-design" begin
        TASOPT.fly_mission!(ac, 2; printTO=false)

        @test ac.parm[imPFEI, 2] ≈ PFEI
        @test ac.parm[imWfuel, 2] ≈ Wfuel
        @test ac.parm[imWTO, 2] ≈ WTO
    end
end

@testset "Default sizing" verbose=true begin
    ac = load_default_model()
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius == 1.9558
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_sized.jl"))
    # Fuselage
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))

    size_aircraft!(ac; printiter=false);

    @testset "Fuselage" begin
        @test  check_struct_equivalence(ac_test.fuselage, ac.fuselage)
    end

    @testset "Wing" begin
        @test  check_struct_equivalence(ac_test.wing, ac.wing)
    end

    @testset "Htail" begin
        @test  check_struct_equivalence(ac_test.htail, ac.htail)
    end

    @testset "Vtail" begin
        @test  check_struct_equivalence(ac_test.vtail, ac.vtail)
    end

    @testset "Geometry" begin
        for i in eachindex(parg)
            @test parg[i] ≈ ac.parg[i] 
        end
    end

    @testset "Aero" begin
        for i in eachindex(para)
            @test para[i] ≈ ac.para[i] 
        end
    end

    @testset "Propulsion (typed EngineState)" begin
        import TOML
        _baseline_path = joinpath(TASOPT.__TASOPTroot__, "../test/fixtures/default_sized_engine_state.toml")
        _bl = TOML.parsefile(_baseline_path)

        # Station field names (mirrors _station_to_dict in engine_harness.jl)
        _st_fields = (:Tt, :ht, :pt, :cpt, :Rt, :Ts, :ps, :cps, :Rs, :u, :A, :mdot)
        # DesignState scalar fields with bare-pare backing
        _ds_fields = (:pi_fan_des, :pi_lpc_des, :pi_hpc_des, :pi_hpt_des, :pi_lpt_des,
                      :mb_fan_des, :mb_lpc_des, :mb_hpc_des, :mb_hpt_des, :mb_lpt_des,
                      :Nb_fan_des, :Nb_lpc_des, :Nb_hpc_des, :Nb_hpt_des, :Nb_lpt_des,
                      :A2, :A25, :A8, :A18,
                      :fc, :ruc, :M4a,
                      :pid, :pib, :pifn, :pitn,
                      :epolf, :epollc, :epolhc, :epolht, :epollt,
                      :pifK, :epfK, :M2, :M25, :epsl, :epsh, :etab,
                      :dTstrk, :Mtexit, :StA, :efilm, :tfilm, :fc0, :dehtdfc)
        # EngineState scalar fields with bare-pare backing
        _eng_fields = (:M0, :T0, :p0, :a0, :rho0, :mu0, :Tfuel, :Tfuel_tank,
                       :RadCoolantT, :RadCoolantP, :Qradiator, :hfuel,
                       :ff, :mofft, :Pofft, :Phiinl, :Kinl,
                       :Nf, :N1, :N2, :Nbf, :Nblc, :Nbhc,
                       :epf, :eplc, :ephc, :epht, :eplt,
                       :TSFC, :Fe, :Fsp, :BPR, :mfuel,
                       :Pfan, :TSEC, :mfan, :Pfanmax,
                       :mbf, :mblc, :mbhc, :pif, :pilc, :pihc,
                       :etaf, :etalc, :etahc, :etaht, :etalt,
                       :eta_thermal, :eta_prop, :eta_overall)
        # Station symbols in TOML_STATION_ORDER order (ARP755 field names)
        _stations = (:st0, :st2, :st12, :st2a, :st2ac, :st13, :st25, :st25c,
                     :st3, :st4, :st4a, :st41, :st45, :st5, :st5c,
                     :st8, :st9, :st18, :st19, :st25off)

        for _bl_pt in _bl["points"]
            _ip  = _bl_pt["ip"]
            _eng = ac.missions[1].points[_ip].engine
            _ds  = _eng.design
            _bl_ds = _bl_pt["design"]

            # Scalar checks
            # rtol=1e-10 atol=1e-14: baseline captured on macOS; CI runs on
            # Linux x86_64 with a different libm/BLAS, and transcendental
            # functions legitimately differ by a few ULPs across platforms.
            # atol backstop handles zero-valued baseline fields (Fsp, Phiinl,
            # Qradiator, ...) that rtol cannot compare meaningfully.
            for _f in _eng_fields
                @test getfield(_eng, _f) ≈ _bl_pt[String(_f)] rtol=1e-10 atol=1e-14
            end
            # DesignState scalar checks
            for _f in _ds_fields
                @test getfield(_ds, _f) ≈ _bl_ds[String(_f)] rtol=1e-10 atol=1e-14
            end
            # DesignState vector fields (epsrow, Tmrow)
            for _i in 1:4
                @test _ds.epsrow[_i] ≈ _bl_ds["epsrow"][_i] rtol=1e-10 atol=1e-14
                @test _ds.Tmrow[_i]  ≈ _bl_ds["Tmrow"][_i]  rtol=1e-10 atol=1e-14
            end
            # FlowStation checks
            _bl_sts = _bl_pt["stations"]
            for _stfld in _stations
                _st    = getfield(_eng, _stfld)
                _bl_st = _bl_sts[String(_stfld)]
                for _sf in _st_fields
                    @test getproperty(_st, _sf) ≈ _bl_st[String(_sf)] rtol=1e-10 atol=1e-14
                end
            end
        end
    end

    test_ac_off_design(ac, 1.0869638391729122, 153128.29535348987,  769359.1150444464)
    
    @test ac.parm[imPFEI] ≈ 0.945758611404728 rtol=1e-4
end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 3.0988
    

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.1903760871373523 rtol=1e-4

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 1.5113

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 0.8483560952994892 rtol=1e-4

end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.0076619899926231 rtol=1e-4

end
