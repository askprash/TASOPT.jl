using StaticArrays

@testset "Ducted fan" begin
    @testset "Thrust from ROC" begin
        ac = TASOPT.load_default_model()
        ip = TASOPT.ipclimb1
        imission = 1

        ac.wing.layout.S = 100.0
        ac.parg[igWMTO] = 1e5
        ac.para[iaCD,ip,imission] = 0.04
        ac.para[iaCL,ip,imission] = 0.5
        ac.para[iafracW,ip,imission] = 1.0
        ac.missions[imission].points[ip].engine.rho0 = 1.2

        ROC = 1.0 #m/s
        
        ac.para[iaROCdes,ip,1] = ROC

        TASOPT.engine.calculate_thrust_from_ROC!(ac, ip, imission)

        @test ac.missions[imission].points[ip].engine.Fe ≈ 4865.490242567016 #Check that the thrust is set to zero when the climb rate is zero
    end
    @testset "Ducted fan models" begin
        #__ Test ducted fan sizing function __
        gee = TASOPT.gee
        M0 = 0.8
        h = 11_000.0 #m
        atmos_state = TASOPT.atmos(h)
        T0 = atmos_state.T
        p0 = atmos_state.p
        ρ0 = atmos_state.ρ
        a0 = atmos_state.a
        μ0 = atmos_state.μ
        M2 = 0.6
        Kinl = 0
        iBLIc = 0
        Phiinl = 0
        pi_fan_des = 1.5
        pid = 1.0
        pifn = 1.0
        epf0 = 0.9
        Δh_radiator = 0
        Δp_radiator = 0
        Fe = 1e4

        out_size = TASOPT.ductedfansize!(gee, M0, T0, p0, a0, M2,
            Fe, Phiinl, Kinl, iBLIc,
            pi_fan_des,
            pid, pifn, 
            epf0,
            Δh_radiator,
            Δp_radiator
            )

        TSEC, Fsp, Pfan, mfan,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, A2,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        epf,
        etaf,
        Lconv = out_size

        out_size_check = (327.11210797234, 0.4556101382778294, 3.2711210797234005e6, 92.71080094936762, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 281.0643310206969, -49272.85600965954, 51990.120356952786, 1005.7936500187604, 287.482334, 281.0643310206969, -49272.85600965954, 51990.120356952786, 1005.7936500187604, 287.482334, 236.74253162629373, 229.40455573526395, 182.38072225868564, 27167.2752649851, 1004.3726660891869, 287.482334, 1.234009545340459, 234.14733320991527, 307.08902200062846, 27455.81621405862, 1004.4531614909928, 287.482334, 0.7401709919109498, 221.89872511133328, 344.82581602886535, 22756.739147503034, 1004.2637990691651, 287.482334, 0.7536791452295546, 0.8693488619036531, 0.8616487982671662, true)
        
        for (i,item) in enumerate(out_size) 
            @test item ≈ out_size_check[i]
        end

        mb_fan_des = mfan * sqrt(Tt2 / TASOPT.Tref) / (pt2 / TASOPT.pref)
        Nf = 1.0 #Arbitrarily set to 1 as only ratios matter
        Nbf = Nf / sqrt(Tt2 / TASOPT.Tref)
        Nb_fan_des = Nbf

        #__ Test ducted fan operation function __
        #First check that it provides the desired values at the design point
        out_opr_des  = TASOPT.ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                    Phiinl, Kinl, iBLIc,
                    pid, pifn, 
                    pi_fan_des, 
                    mb_fan_des,
                    Nb_fan_des, 
                    A2, A7,
                    epf0,
                    Fe, 0,
                    M2, pi_fan_des, 0, 
                    Δh_radiator, Δp_radiator,
                    false)

        Fe2 = out_opr_des[3]
        Pf2 = out_opr_des[4]
        mfan2 = out_opr_des[5]
        pif2 = out_opr_des[6]

        #Test that the operation code returns the design variables
        @test Fe2 ≈ Fe
        @test Pf2 ≈ Pfan
        @test mfan2 ≈ mfan
        @test pif2 ≈ pi_fan_des

        #Now check conditions at takeoff 
        atmos_state = TASOPT.atmos(0.0)
        T0 = atmos_state.T
        p0 = atmos_state.p
        ρ0 = atmos_state.ρ
        a0 = atmos_state.a
        μ0 = atmos_state.μ
        M0 = 0
        Pf_takeoff = 8e6

        out_opr_takeoff  = TASOPT.ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                    Phiinl, Kinl, iBLIc,
                    pid, pifn, 
                    pi_fan_des, 
                    mb_fan_des,
                    Nb_fan_des, 
                    A2, A7,
                    epf0,
                    0, Pf_takeoff,
                    M2, pi_fan_des, 0, 
                    Δh_radiator, Δp_radiator,
                    true)
        
        out_opr_check = (140.74270917523936, 0.0, 56841.31026666052, 8.000000000000149e6, 225.39505852577363, 1.4330295431713083, 225.39505852577363, 1.0012505100170295, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 323.4458621045641, -6600.376081204027, 145194.55331411696, 1007.9218437076298, 287.482334, 323.4458621045641, -6600.376081204027, 145194.55331411696, 1007.9218437076298, 287.482334, 0.0, 273.9002717202924, 169.61337292744685, 84794.83629183592, 1005.535918910043, 287.482334, 0.5107847532243659, 291.8712909624618, 252.18525480744196, 101320.0, 1006.2219870034098, 287.482334, 0.7357958274607429, 291.8712909624618, 252.18525480744196, 101320.0, 1006.2219870034098, 287.482334, 0.7357958274607429, 0.7401709919114717, 0.8903649191567354, 0.8846551852455233)

        for (i,item) in enumerate(out_opr_takeoff) 
            @test item ≈ out_opr_check[i]
        end
    end

    @testset "Ducted fan with fuel cell" begin
        ac = TASOPT.load_default_model()
        parg = ac.parg

        #First, test the fuel cell + ducted fan design
        #Create a ducted fan with fuel cell model
        modelname = "fuel_cell_with_ducted_fan"
        engineweightname = "nasa"

        enginecalc! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
        engineweight! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
        enginemodel = TASOPT.engine.FuelCellDuctedFan(modelname, enginecalc!, engineweightname, engineweight!, false)

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
        fcdata.thickness_anode  = 250e-6
        fcdata.thickness_cathode  = 250e-6
        fcdata.design_voltage = 200.0

        ac.para[iaROCdes, ipclimb1:ipclimbn,:] .= 500 * ft_to_m / 60
        engdata = fcdata

        engine = TASOPT.engine.Engine(enginemodel, engdata, Vector{TASOPT.engine.HeatExchanger}())

        ac.engine = engine

        # Seed per-point typed state: Pfanmax and radiator fields.
        for im in 1:length(ac.missions), ip in 1:length(ac.missions[im].points)
            eng_ip = ac.missions[im].points[ip].engine
            eng_ip.Pfanmax      = 10e6
            eng_ip.RadCoolantT  = engine.data.FC_temperature[ip, im]
            eng_ip.RadCoolantP  = engine.data.FC_pressure[ip, im]
            eng_ip.Qradiator    = engine.data.FC_heat[ip, im]
        end

        # Design-point inputs at ipcruise1 written directly to typed state.
        let eng = ac.missions[1].points[ipcruise1].engine
            eng.Fe             = 16981.808185580507
            eng.Fsp            = 0.5268888878557653
            eng.pif            = 1.685
            eng.design.pid     = 0.998
            eng.design.pifn    = 0.98
            eng.design.epolf   = 0.8948
            eng.design.pifK    = 1.685
            eng.design.epfK    = -0.077
            eng.design.M2      = 0.6
            eng.M0             = 0.8
            eng.st0.Tt         = 247.50937622665612
            eng.st0.ht         = -83004.69052610746
            eng.st0.pt         = 36436.066979351635
            eng.st0.cpt        = 1004.734911318952
            eng.st0.Rt         = 287.482334
            eng.p0             = 23922.608843328788
            eng.a0             = 296.8557888469756
            eng.rho0           = 0.3800541947033053
            eng.mu0            = 1.4294279408521106e-5
            eng.T0             = 219.43067572699252
            eng.st0.u          = 237.4846310775805
        end
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
        # FC design uses ip_fcdes=iprotate for Pfanmax (already set above in the loop).
        engine.enginecalc!(ac, "design", 1, ipcruise1, true)

        @test ac.missions[1].points[ipcruise1].engine.mfuel ≈ 0.0005481461619779067
        @test ac.missions[1].points[ipcruise1].engine.Pfan  ≈ 6.043228014711607e6
        @test ac.engine.data.number_cells ≈ 265.5192500533146
        @test ac.engine.data.area_cell ≈ 5.0
        @test ac.engine.data.FC_heat[ipcruise1,1] ≈ 2.6917275229945835e6

        #Next, test the off-design performance — write off-design inputs directly to typed state.
        let eng = ac.missions[1].points[iprotate].engine
            eng.Fe           = 20e3
            eng.Fsp          = 3.25725288027002
            eng.pif          = 1.7
            eng.design.pid   = 0.998
            eng.design.pifn  = 0.98
            eng.design.epolf = 0.8948
            eng.M0           = 0.21721987853710417
            eng.st0.Tt       = 290.9134582811252
            eng.st0.ht       = -39363.02002852104
            eng.st0.pt       = 104698.14961489625
            eng.st0.cpt      = 1006.18147266869
            eng.st0.Rt       = 287.482334
            eng.p0           = 101320.0
            eng.a0           = 340.2074661144284
            eng.rho0         = 1.2255627040761317
            eng.mu0          = 1.78e-5
            eng.T0           = 288.2
            eng.st0.u        = 73.89982446679213
        end
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
        engine.enginecalc!(ac, "off_design", 1, iprotate, true)

        @test ac.missions[1].points[iprotate].engine.mfuel ≈ 0.0010447183729991173
        @test ac.missions[1].points[iprotate].engine.Pfan  ≈ 1.0000000000384096e7
        @test ac.engine.data.FC_heat[iprotate,1] ≈ 6.64805697887804e6
    end

    # =========================================================================
    # run_ducted_fan_design_point / pare_to_ducted_fan_state! (tasopt-3gk)
    # =========================================================================
    @testset "Ducted fan harness" begin
        import TOML
        # ------------------------------------------------------------------
        # Recreate the fuel-cell ducted-fan aircraft from scratch.
        # Mirrors the "Ducted fan with fuel cell" setup above but also sets
        # para[iaalt/iaMach] and parm[imaltTO/imT0TO] required by the harness.
        # ------------------------------------------------------------------
        ac_h = TASOPT.load_default_model()
        parg_h = ac_h.parg

        enginecalc_h! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
        engineweight_h! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
        enginemodel_h = TASOPT.engine.FuelCellDuctedFan(
            "fuel_cell_with_ducted_fan", enginecalc_h!, "nasa", engineweight_h!, false)

        fcdata_h = TASOPT.engine.FuelCellDuctedFanData(2)
        fcdata_h.type = "HT-PEMFC"
        fcdata_h.current_density[iprotate,:] .= 1e4
        fcdata_h.FC_temperature .= 453.15
        fcdata_h.FC_pressure .= 3e5
        fcdata_h.water_concentration_anode .= 0.1
        fcdata_h.water_concentration_cathode .= 0.1
        fcdata_h.λ_H2 .= 3.0
        fcdata_h.λ_O2 .= 3.0
        fcdata_h.thickness_membrane = 100e-6
        fcdata_h.thickness_anode    = 250e-6
        fcdata_h.thickness_cathode  = 250e-6
        fcdata_h.design_voltage = 200.0
        ac_h.para[iaROCdes, ipclimb1:ipclimbn, :] .= 500 * ft_to_m / 60

        engine_h = TASOPT.engine.Engine(enginemodel_h, fcdata_h,
                                        Vector{TASOPT.engine.HeatExchanger}())
        ac_h.engine = engine_h

        # Seed per-point typed state: Pfanmax and radiator fields.
        for im in 1:length(ac_h.missions), ip in 1:length(ac_h.missions[im].points)
            eng_ip = ac_h.missions[im].points[ip].engine
            eng_ip.Pfanmax     = 10e6
            eng_ip.RadCoolantT = engine_h.data.FC_temperature[ip, im]
            eng_ip.RadCoolantP = engine_h.data.FC_pressure[ip, im]
            eng_ip.Qradiator   = engine_h.data.FC_heat[ip, im]
        end

        # Design-point engine inputs at ipcruise1 written directly to typed state.
        let eng = ac_h.missions[1].points[ipcruise1].engine
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
        parg_h[igneng] = 2
        ac_h.wing.layout.S = 81.25043040696103

        # Fields required by run_ducted_fan_design_point: altitude and Mach.
        # 10668 m ≈ FL350 standard-day matches the manually-set T0/p0 above.
        ac_h.para[iaalt,  ipcruise1, 1] = 10668.0
        ac_h.para[iaMach, ipcruise1, 1] = 0.8
        # parm defaults from default_input.toml: imaltTO=0 (sea level), imT0TO=288.2K
        # FC design uses ip_fcdes=iprotate for Pfanmax (already seeded above in the loop).

        # ------------------------------------------------------------------
        # run_ducted_fan_design_point
        # ------------------------------------------------------------------
        df = TASOPT.engine.run_ducted_fan_design_point(ac_h)

        # Ambient scalars: physically plausible cruise conditions
        @test df.M0 ≈ 0.8       rtol=1e-8
        @test df.T0 > 200.0                 # stratosphere static temp, K
        @test df.p0 > 1e4                   # cruise static pressure, Pa
        @test df.a0 > 280.0                 # speed of sound, m/s

        # Performance: sensible values at design thrust
        @test df.Fe   ≈ 16981.808185580507 rtol=1e-10
        @test df.Pfan > 0.0
        @test df.mfan > 0.0
        @test 0.0 < df.etaf < 1.0

        # Fan pressure ratio > 1 (compression work)
        @test df.pif > 1.0

        # Design areas are positive
        @test df.A2 > 0.0
        @test df.A7 > 0.0

        # Fan mass flow stored in st2.mdot
        @test df.st2.mdot > 0.0

        # Thermodynamic invariant: fan adds enthalpy → Tt rises
        @test df.st21.Tt > df.st2.Tt    # fan exit hotter than inlet

        # Fan nozzle is adiabatic: Tt7 ≈ Tt21
        @test df.st7.Tt ≈ df.st21.Tt rtol=1e-8

        # Nozzle ideally expanded at design point: p_exit ≈ p0
        @test df.st8.ps ≈ df.p0 rtol=1e-6

        # Round-trip consistency: DuctedFanState values match typed EngineState
        let eng = ac_h.missions[1].points[ipcruise1].engine
            @test df.Fe     ≈ 16981.808185580507 rtol=1e-12   # Fe: test-set input
            @test df.Pfan   ≈ eng.Pfan          rtol=1e-12
            @test df.A2     ≈ eng.design.A2     rtol=1e-12
            @test df.A7     ≈ eng.design.A7     rtol=1e-12
            @test df.pif    ≈ 1.685             rtol=1e-12   # pif: test-set input
            @test df.etaf   ≈ eng.etaf          rtol=1e-12
            @test df.st2.Tt ≈ eng.st2.Tt       rtol=1e-12
        end

        # ------------------------------------------------------------------
        # pare_to_ducted_fan_state! round-trip: fresh state matches harness
        # ------------------------------------------------------------------
        df2 = TASOPT.engine.DuctedFanState{Float64}()
        TASOPT.engine.pare_to_ducted_fan_state!(df2, ac_h.missions[1].points[ipcruise1].engine)
        @test df2.Fe   ≈ df.Fe   rtol=1e-12
        @test df2.pif  ≈ df.pif  rtol=1e-12
        @test df2.A2   ≈ df.A2   rtol=1e-12
        @test df2.st21.Tt ≈ df.st21.Tt rtol=1e-12

        # ------------------------------------------------------------------
        # run_ducted_fan_sweep — single-point off-design at design conditions.
        # At the design point, off-design should reproduce the design-point state.
        # ------------------------------------------------------------------
        states = TASOPT.engine.run_ducted_fan_sweep(ac_h; ip_range=ipcruise1:ipcruise1)

        @test haskey(states, ipcruise1)
        df_sweep = states[ipcruise1]
        @test df_sweep.Fe   ≈ df.Fe   rtol=1e-6
        @test df_sweep.pif  ≈ df.pif  rtol=1e-6
        @test df_sweep.etaf ≈ df.etaf rtol=1e-6

        # ------------------------------------------------------------------
        # write_ducted_fan_sweep_toml — valid TOML with expected structure.
        # ------------------------------------------------------------------
        buf = IOBuffer()
        TASOPT.engine.write_ducted_fan_sweep_toml(buf, states, ipcruise1:ipcruise1, ac_h)
        toml_str = String(take!(buf))
        toml_data = TOML.parse(toml_str)

        @test haskey(toml_data, "points")
        @test length(toml_data["points"]) == 1
        pt = toml_data["points"][1]
        @test pt["ip"] == ipcruise1
        @test pt["Fe_N"] ≈ df.Fe rtol=1e-10

        # Station subtables present
        @test haskey(pt, "stations")
        @test haskey(pt["stations"], "st2")
        @test haskey(pt["stations"], "st7")
        @test pt["stations"]["st2"]["Tt"] ≈ df.st2.Tt rtol=1e-10
        @test pt["stations"]["st7"]["Tt"] ≈ df.st7.Tt rtol=1e-10

        # ------------------------------------------------------------------
        # Regression baseline: compare against pinned fixture.
        # If this fails, regenerate with:
        #   julia --project=. test/generate_ducted_fan_baseline.jl
        # ------------------------------------------------------------------
        baseline_path = TASOPT.engine.DUCTED_FAN_BASELINE_PATH
        if isfile(baseline_path)
            bl = TOML.parsefile(baseline_path)
            bl_pt = bl["points"][1]

            @test bl_pt["ip"] == ipcruise1
            @test bl_pt["Fe_N"]        ≈ df.Fe   rtol=1e-10
            @test bl_pt["etaf"]        ≈ df.etaf rtol=1e-10
            @test bl_pt["Pfan_W"]      ≈ df.Pfan rtol=1e-10
            @test bl_pt["pif"]         ≈ df.pif  rtol=1e-10
            @test bl_pt["mfan_kg_s"]   ≈ df.mfan rtol=1e-10

            bl_st2 = bl_pt["stations"]["st2"]
            @test bl_st2["Tt"] ≈ df.st2.Tt rtol=1e-10
            @test bl_st2["pt"] ≈ df.st2.pt rtol=1e-10

            bl_st7 = bl_pt["stations"]["st7"]
            @test bl_st7["Tt"] ≈ df.st7.Tt rtol=1e-10
            @test bl_st7["u"]  ≈ df.st7.u  rtol=1e-10
        end
    end  # Ducted fan harness
end
