using StaticArrays

# Test-local shim: test fixtures seed typed state via bare pare, so we need to
# sync pare → typed EngineState before each enginecalc! call that relies on it.
# Production callers populate typed state directly; the src-side pre-sync was
# removed in tasopt-j9l.45.14.6.5.
# pare_to_engine_state! was deleted in tasopt-j9l.45.14.6; body inlined here.
function _sync_pare_to_engine_state!(ac, imission, ip)
    eng  = ac.missions[imission].points[ip].engine
    pare = view(ac.pare, :, ip, imission)
    eng.M0    = pare[ieM0];   eng.T0    = pare[ieT0];   eng.p0    = pare[iep0]
    eng.a0    = pare[iea0];   eng.rho0  = pare[ierho0]; eng.mu0   = pare[iemu0]
    eng.Tfuel      = pare[ieTfuel];   eng.Tfuel_tank = pare[ieTft]
    eng.RadCoolantT = pare[ieRadiatorCoolantT]; eng.RadCoolantP = pare[ieRadiatorCoolantP]
    eng.Qradiator   = pare[ieRadiatorHeat];     eng.hfuel  = pare[iehfuel]
    eng.ff = pare[ieff]; eng.mofft = pare[iemofft]; eng.Pofft = pare[iePofft]
    eng.Phiinl = pare[iePhiinl]; eng.Kinl = pare[ieKinl]
    eng.Nf = pare[ieNf]; eng.N1 = pare[ieN1]; eng.N2 = pare[ieN2]
    eng.Nbf = pare[ieNbf]; eng.Nblc = pare[ieNblc]; eng.Nbhc = pare[ieNbhc]
    eng.epf = pare[ieepf]; eng.eplc = pare[ieeplc]; eng.ephc = pare[ieephc]
    eng.epht = pare[ieepht]; eng.eplt = pare[ieeplt]
    # Station 0
    eng.st0.Tt = pare[ieTt0]; eng.st0.ht = pare[ieht0]; eng.st0.pt = pare[iept0]
    eng.st0.cpt = pare[iecpt0]; eng.st0.Rt = pare[ieRt0]; eng.st0.u = pare[ieu0]
    # Station 18
    eng.st18.Tt = pare[ieTt18]; eng.st18.ht = pare[ieht18]; eng.st18.pt = pare[iept18]
    eng.st18.cpt = pare[iecpt18]; eng.st18.Rt = pare[ieRt18]
    # Station 19
    eng.st19.Tt = pare[ieTt19]; eng.st19.ht = pare[ieht19]; eng.st19.pt = pare[iept19]
    eng.st19.cpt = pare[iecpt19]; eng.st19.Rt = pare[ieRt19]
    # Station 2
    eng.st2.Tt = pare[ieTt2]; eng.st2.ht = pare[ieht2]; eng.st2.pt = pare[iept2]
    eng.st2.cpt = pare[iecpt2]; eng.st2.Rt = pare[ieRt2]
    eng.st2.ps = pare[iep2]; eng.st2.Ts = pare[ieT2]; eng.st2.Rs = pare[ieR2]
    eng.st2.cps = pare[iecp2]; eng.st2.u = pare[ieu2]
    eng.st2.A = pare[ieA2]; eng.st2.mdot = pare[iemcore]
    # Station 21
    eng.st21.Tt = pare[ieTt21]; eng.st21.ht = pare[ieht21]; eng.st21.pt = pare[iept21]
    eng.st21.cpt = pare[iecpt21]; eng.st21.Rt = pare[ieRt21]
    # Station 25
    eng.st25.Tt = pare[ieTt25]; eng.st25.ht = pare[ieht25]; eng.st25.pt = pare[iept25]
    eng.st25.cpt = pare[iecpt25]; eng.st25.Rt = pare[ieRt25]
    eng.st25.ps = pare[iep25]; eng.st25.Ts = pare[ieT25]; eng.st25.Rs = pare[ieR25]
    eng.st25.cps = pare[iecp25]; eng.st25.u = pare[ieu25]; eng.st25.A = pare[ieA25]
    # Station 3
    eng.st3.Tt = pare[ieTt3]; eng.st3.ht = pare[ieht3]; eng.st3.pt = pare[iept3]
    eng.st3.cpt = pare[iecpt3]; eng.st3.Rt = pare[ieRt3]
    # Station 4
    eng.st4.Tt = pare[ieTt4]; eng.st4.ht = pare[ieht4]; eng.st4.pt = pare[iept4]
    eng.st4.cpt = pare[iecpt4]; eng.st4.Rt = pare[ieRt4]
    # Station 41
    eng.st41.Tt = pare[ieTt41]; eng.st41.ht = pare[ieht41]; eng.st41.pt = pare[iept41]
    eng.st41.cpt = pare[iecpt41]; eng.st41.Rt = pare[ieRt41]
    # Station 45
    eng.st45.Tt = pare[ieTt45]; eng.st45.ht = pare[ieht45]; eng.st45.pt = pare[iept45]
    eng.st45.cpt = pare[iecpt45]; eng.st45.Rt = pare[ieRt45]
    # Station 49
    eng.st49.Tt = pare[ieTt49]; eng.st49.ht = pare[ieht49]; eng.st49.pt = pare[iept49]
    eng.st49.cpt = pare[iecpt49]; eng.st49.Rt = pare[ieRt49]
    # Station 5
    eng.st5.Tt = pare[ieTt5]; eng.st5.ht = pare[ieht5]; eng.st5.pt = pare[iept5]
    eng.st5.cpt = pare[iecpt5]; eng.st5.Rt = pare[ieRt5]
    eng.st5.ps = pare[iep5]; eng.st5.Ts = pare[ieT5]; eng.st5.Rs = pare[ieR5]
    eng.st5.cps = pare[iecp5]; eng.st5.u = pare[ieu5]; eng.st5.A = pare[ieA5]
    # Station 6
    eng.st6.ps = pare[iep6]; eng.st6.Ts = pare[ieT6]; eng.st6.Rs = pare[ieR6]
    eng.st6.cps = pare[iecp6]; eng.st6.u = pare[ieu6]; eng.st6.A = pare[ieA6]
    # Station 7
    eng.st7.Tt = pare[ieTt7]; eng.st7.ht = pare[ieht7]; eng.st7.pt = pare[iept7]
    eng.st7.cpt = pare[iecpt7]; eng.st7.Rt = pare[ieRt7]
    eng.st7.ps = pare[iep7]; eng.st7.Ts = pare[ieT7]; eng.st7.Rs = pare[ieR7]
    eng.st7.cps = pare[iecp7]; eng.st7.u = pare[ieu7]; eng.st7.A = pare[ieA7]
    # Station 8
    eng.st8.ps = pare[iep8]; eng.st8.Ts = pare[ieT8]; eng.st8.Rs = pare[ieR8]
    eng.st8.cps = pare[iecp8]; eng.st8.u = pare[ieu8]; eng.st8.A = pare[ieA8]
    # Station 9
    eng.st9.Tt = pare[ieTt9]; eng.st9.pt = pare[iept9]; eng.st9.u = pare[ieu9]
    eng.st9.A  = pare[ieA9]
    # Cooling design state
    eng.design.epsrow = SVector{4,Float64}(pare[ieepsc1], pare[ieepsc2], pare[ieepsc3], pare[ieepsc4])
    eng.design.Tmrow  = SVector{4,Float64}(pare[ieTmet1], pare[ieTmet2], pare[ieTmet3], pare[ieTmet4])
    eng.design.fc  = pare[iefc]; eng.design.ruc = pare[ieruc]; eng.design.M4a = pare[ieM4a]
    # Design map anchors
    eng.design.pifD = pare[iepifD]; eng.design.pilcD = pare[iepilcD]
    eng.design.pihcD = pare[iepihcD]; eng.design.pihtD = pare[iepihtD]; eng.design.piltD = pare[iepiltD]
    eng.design.mbfD = pare[iembfD]; eng.design.mblcD = pare[iemblcD]
    eng.design.mbhcD = pare[iembhcD]; eng.design.mbhtD = pare[iembhtD]; eng.design.mbltD = pare[iembltD]
    eng.design.NbfD = pare[ieNbfD]; eng.design.NblcD = pare[ieNblcD]
    eng.design.NbhcD = pare[ieNbhcD]; eng.design.NbhtD = pare[ieNbhtD]; eng.design.NbltD = pare[ieNbltD]
    eng.design.A2 = pare[ieA2]; eng.design.A25 = pare[ieA25]
    eng.design.A5 = pare[ieA5]; eng.design.A7  = pare[ieA7]
    # Performance rollup
    eng.TSFC = pare[ieTSFC]; eng.Fe = pare[ieFe]; eng.Fsp = pare[ieFsp]
    eng.BPR  = pare[ieBPR];  eng.mfuel = pare[iemfuel]
    # Ducted-fan outputs
    eng.Pfan = pare[iePfan]; eng.TSEC = pare[ieTSEC]; eng.mfan = pare[iemfan]
    eng.Pfanmax = pare[iePfanmax]
    # Compressor map operating points
    eng.mbf = pare[iembf]; eng.mblc = pare[iemblc]; eng.mbhc = pare[iembhc]
    eng.pif = pare[iepif]; eng.pilc = pare[iepilc]; eng.pihc = pare[iepihc]
    # Adiabatic efficiencies
    eng.etaf = pare[ieetaf]; eng.etalc = pare[ieetalc]; eng.etahc = pare[ieetahc]
    eng.etaht = pare[ieetaht]; eng.etalt = pare[ieetalt]
    # Design-constant scalars
    eng.design.pid = pare[iepid]; eng.design.pib = pare[iepib]
    eng.design.pifn = pare[iepifn]; eng.design.pitn = pare[iepitn]
    eng.design.epolf = pare[ieepolf]; eng.design.epollc = pare[ieepollc]
    eng.design.epolhc = pare[ieepolhc]; eng.design.epolht = pare[ieepolht]; eng.design.epollt = pare[ieepollt]
    eng.design.pifK = pare[iepifK]; eng.design.epfK = pare[ieepfK]
    eng.design.M2 = pare[ieM2]; eng.design.M25 = pare[ieM25]
    eng.design.epsl = pare[ieepsl]; eng.design.epsh = pare[ieepsh]; eng.design.etab = pare[ieetab]
    eng.design.dTstrk = pare[iedTstrk]; eng.design.Mtexit = pare[ieMtexit]; eng.design.StA = pare[ieStA]
    eng.design.efilm = pare[ieefilm]; eng.design.tfilm = pare[ietfilm]
    eng.design.fc0 = pare[iefc0]; eng.design.dehtdfc = pare[iedehtdfc]
    return eng
end

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
        ac.pare[ierho0,ip,imission] = 1.2
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
        pifD = 1.5
        pid = 1.0
        pifn = 1.0
        epf0 = 0.9
        Δh_radiator = 0
        Δp_radiator = 0
        Fe = 1e4

        out_size = TASOPT.ductedfansize!(gee, M0, T0, p0, a0, M2,
            Fe, Phiinl, Kinl, iBLIc,
            pifD,
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

        mbfD = mfan * sqrt(Tt2 / TASOPT.Tref) / (pt2 / TASOPT.pref)
        Nf = 1.0 #Arbitrarily set to 1 as only ratios matter
        Nbf = Nf / sqrt(Tt2 / TASOPT.Tref)
        NbfD = Nbf

        #__ Test ducted fan operation function __
        #First check that it provides the desired values at the design point
        out_opr_des  = TASOPT.ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                    Phiinl, Kinl, iBLIc,
                    pid, pifn, 
                    pifD, 
                    mbfD,
                    NbfD, 
                    A2, A7,
                    epf0,
                    Fe, 0,
                    M2, pifD, 0, 
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
        @test pif2 ≈ pifD

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
                    pifD, 
                    mbfD,
                    NbfD, 
                    A2, A7,
                    epf0,
                    0, Pf_takeoff,
                    M2, pifD, 0, 
                    Δh_radiator, Δp_radiator,
                    true)
        
        out_opr_check = (140.74270917523936, 0.0, 56841.31026666052, 8.000000000000149e6, 225.39505852577363, 1.4330295431713083, 225.39505852577363, 1.0012505100170295, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 323.4458621045641, -6600.376081204027, 145194.55331411696, 1007.9218437076298, 287.482334, 323.4458621045641, -6600.376081204027, 145194.55331411696, 1007.9218437076298, 287.482334, 0.0, 273.9002717202924, 169.61337292744685, 84794.83629183592, 1005.535918910043, 287.482334, 0.5107847532243659, 291.8712909624618, 252.18525480744196, 101320.0, 1006.2219870034098, 287.482334, 0.7357958274607429, 291.8712909624618, 252.18525480744196, 101320.0, 1006.2219870034098, 287.482334, 0.7357958274607429, 0.7401709919114717, 0.8903649191567354, 0.8846551852455233)

        for (i,item) in enumerate(out_opr_takeoff) 
            @test item ≈ out_opr_check[i]
        end
    end

    @testset "Ducted fan with fuel cell" begin
        ac = TASOPT.load_default_model()
        pare = ac.pare
        parg = ac.parg

        #First, test the fuel cell + ducted fan design
        #Create a ducted fan with fuel cell model
        modelname = "fuel_cell_with_ducted_fan"
        engineweightname = "nasa"

        enginecalc! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
        engineweight! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
        enginemodel = TASOPT.engine.FuelCellDuctedFan(modelname, enginecalc!, engineweightname, engineweight!, false)
        pare[iePfanmax,:,:] .= 10e6

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
        pare[ieRadiatorepsilon,:,:] .= 0.7
        pare[ieRadiatorMp,:,:] .= 0.12
        pare[ieDi,:,:] .= 0.4

        para[iaROCdes, ipclimb1:ipclimbn,:] .= 500 * ft_to_m / 60
        engdata = fcdata

        engine = TASOPT.engine.Engine(enginemodel, engdata, Vector{TASOPT.engine.HeatExchanger}())

        ac.engine = engine

        #Prepare the pare object
        pare[ieRadiatorCoolantT,:,:] = engine.data.FC_temperature[:,:]
        pare[ieRadiatorCoolantP,:,:] = engine.data.FC_pressure[:,:]
        pare[ieRadiatorHeat,:,:] = engine.data.FC_heat[:,:]

        pare[ieFe,ipcruise1,1] = 16981.808185580507
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
    
        pare[ieFsp,ipcruise1,1] = 0.5268888878557653
        pare[iepif,ipcruise1,1] = 1.685
        pare[iepid,ipcruise1,1] = 0.998
        pare[iepifn,ipcruise1,1] = 0.98
        
        pare[ieepolf,ipcruise1,1] = 0.8948
        
        pare[iepifK,ipcruise1,1] = 1.685
        pare[ieepfK,ipcruise1,1] = -0.077
        pare[ieM2,ipcruise1,1] = 0.6
        pare[ieM0,ipcruise1,1] = 0.8
        pare[ieTt0,ipcruise1,1] = 247.50937622665612
        pare[ieht0,ipcruise1,1] = -83004.69052610746
        pare[iept0,ipcruise1,1] = 36436.066979351635
        pare[iecpt0,ipcruise1,1] = 1004.734911318952
        pare[ieRt0,ipcruise1,1] = 287.482334
        pare[iep0,ipcruise1,1] = 23922.608843328788
        pare[iea0,ipcruise1,1] = 296.8557888469756
        pare[ierho0,ipcruise1,1] = 0.3800541947033053
        pare[iemu0,ipcruise1,1] = 1.4294279408521106e-5
        pare[ieT0,ipcruise1,1] = 219.43067572699252
        pare[ieu0,ipcruise1,1] = 237.4846310775805

        _sync_pare_to_engine_state!(ac, 1, ipcruise1)
        engine.enginecalc!(ac, "design", 1, ipcruise1, true)

        @test ac.missions[1].points[ipcruise1].engine.mfuel ≈ 0.0005481461619779067
        @test ac.missions[1].points[ipcruise1].engine.Pfan  ≈ 6.043228014711607e6
        @test ac.engine.data.number_cells ≈ 265.5192500533146
        @test ac.engine.data.area_cell ≈ 5.0
        @test ac.engine.data.FC_heat[ipcruise1,1] ≈ 2.6917275229945835e6

        #Next, test the off-design performance

        pare[ieFe,iprotate,1] = 20e3
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
    
        pare[ieFsp,iprotate,1] = 3.25725288027002
        pare[iepif,iprotate,1] = 1.718876825147152
        pare[iepid,iprotate,1] = 0.998
        pare[iepifn,iprotate,1] = 0.98
        
        pare[ieepolf,iprotate,1] = 0.8948
        
        pare[iepif,iprotate,1] = 1.7
        pare[ieM2,iprotate,1] = 0.6
        pare[ieM0,iprotate,1] = 0.21721987853710417
        pare[ieTt0,iprotate,1] = 290.9134582811252
        pare[ieht0,iprotate,1] = -39363.02002852104
        pare[iept0,iprotate,1] = 104698.14961489625
        pare[iecpt0,iprotate,1] = 1006.18147266869
        pare[ieRt0,iprotate,1] = 287.482334
        pare[iep0,iprotate,1] = 101320.0
        pare[iea0,iprotate,1] = 340.2074661144284
        pare[ierho0,iprotate,1] = 1.2255627040761317
        pare[iemu0,iprotate,1] = 1.78e-5
        pare[ieT0,iprotate,1] = 288.2
        pare[ieu0,iprotate,1] = 73.89982446679213
        _sync_pare_to_engine_state!(ac, 1, iprotate)
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
        pare_h = ac_h.pare
        parg_h = ac_h.parg

        enginecalc_h! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
        engineweight_h! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
        enginemodel_h = TASOPT.engine.FuelCellDuctedFan(
            "fuel_cell_with_ducted_fan", enginecalc_h!, "nasa", engineweight_h!, false)
        pare_h[iePfanmax,:,:] .= 10e6

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
        pare_h[ieRadiatorepsilon,:,:] .= 0.7
        pare_h[ieRadiatorMp,:,:]      .= 0.12
        pare_h[ieDi,:,:]              .= 0.4
        ac_h.para[iaROCdes, ipclimb1:ipclimbn, :] .= 500 * ft_to_m / 60

        engine_h = TASOPT.engine.Engine(enginemodel_h, fcdata_h,
                                        Vector{TASOPT.engine.HeatExchanger}())
        ac_h.engine = engine_h

        pare_h[ieRadiatorCoolantT,:,:] = engine_h.data.FC_temperature[:,:]
        pare_h[ieRadiatorCoolantP,:,:] = engine_h.data.FC_pressure[:,:]
        pare_h[ieRadiatorHeat,:,:]     = engine_h.data.FC_heat[:,:]

        # Design-point engine inputs at ipcruise1
        pare_h[ieFe,    ipcruise1, 1] = 16981.808185580507
        parg_h[igneng]  = 2
        ac_h.wing.layout.S = 81.25043040696103
        pare_h[ieFsp,   ipcruise1, 1] = 0.5268888878557653
        pare_h[iepif,   ipcruise1, 1] = 1.685
        pare_h[iepid,   ipcruise1, 1] = 0.998
        pare_h[iepifn,  ipcruise1, 1] = 0.98
        pare_h[ieepolf, ipcruise1, 1] = 0.8948
        pare_h[iepifK,  ipcruise1, 1] = 1.685
        pare_h[ieepfK,  ipcruise1, 1] = -0.077
        pare_h[ieM2,    ipcruise1, 1] = 0.6

        # Fields required by run_ducted_fan_design_point: altitude and Mach.
        # 10668 m ≈ FL350 standard-day matches the manually-set T0/p0 above.
        ac_h.para[iaalt,  ipcruise1, 1] = 10668.0
        ac_h.para[iaMach, ipcruise1, 1] = 0.8
        # parm defaults from default_input.toml: imaltTO=0 (sea level), imT0TO=288.2K

        # Sync bare-pare design inputs → typed EngineState before the harness call
        # (the src-side pre-sync in calculate_fuel_cell_with_ducted_fan! was removed
        # in tasopt-j9l.45.14.6.5; this shim is the test-local replacement).
        _sync_pare_to_engine_state!(ac_h, 1, ipcruise1)

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

        # Round-trip consistency: DuctedFanState values match pare directly
        @test df.Fe   ≈ pare_h[ieFe,   ipcruise1, 1] rtol=1e-12
        @test df.Pfan ≈ pare_h[iePfan, ipcruise1, 1] rtol=1e-12
        @test df.A2   ≈ pare_h[ieA2,   ipcruise1, 1] rtol=1e-12
        @test df.A7   ≈ pare_h[ieA7,   ipcruise1, 1] rtol=1e-12
        @test df.pif  ≈ pare_h[iepif,  ipcruise1, 1] rtol=1e-12
        @test df.etaf ≈ pare_h[ieetaf, ipcruise1, 1] rtol=1e-12
        @test df.st2.Tt ≈ pare_h[ieTt2, ipcruise1, 1] rtol=1e-12

        # ------------------------------------------------------------------
        # pare_to_ducted_fan_state! round-trip: fresh state matches harness
        # ------------------------------------------------------------------
        df2 = TASOPT.engine.DuctedFanState{Float64}()
        TASOPT.engine.pare_to_ducted_fan_state!(df2, view(ac_h.pare, :, ipcruise1, 1))
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
