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

        ROC = 1.0 #m/s
        
        ac.para[iaROCdes,ip,1] = ROC

        TASOPT.engine.calculate_thrust_from_ROC!(ac, ip, imission)

        @test ac.pare[ieFe,ip,1] ≈ 4865.490242567016 #Check that the thrust is set to zero when the climb rate is zero
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

        out_size_check = (327.104580011918, 0.45541741829864213, 3.27104580011918e6, 92.75003358060107, 245.960634581053, -84685.92659213145, 34662.4644792426, 1004.7049518646987, 287.38522541059365, 245.960634581053, -84685.92659213145, 34662.4644792426, 1004.7049518646987, 287.38522541059365, 245.960634581053, -84685.92659213145, 34662.4644792426, 1004.7049518646987, 287.38522541059365, 281.04684612435744, -49418.59919783357, 51993.696718863896, 1005.717678561198, 287.38522541059365, 281.04684612435744, -49418.59919783357, 51993.696718863896, 1005.717678561198, 287.38522541059365, 236.74253162629373, 229.41206161826514, 182.33937393325425, 27172.138088350323, 1004.4049989256071, 287.38522541059365, 1.2342140558319696, 234.15218830271218, 307.0178018111572, 27460.168034167855, 1004.4816061874145, 287.38522541059365, 0.7403038111919118, 221.89725254126807, 344.7858713169373, 22756.739147503034, 1004.2975256478297, 287.38522541059365, 0.7538257097935404, 0.8693488619036531, 0.8616537588268466, true)
        
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
        
        out_opr_check = (140.71891697457957, 0.0, 56850.92077169262, 8.000000000000149e6, 225.47249144652213, 1.4330945576783753, 225.47249144652213, 1.0013069808658677, 288.2, -42223.5408832171, 101320.0, 1006.0037676032673, 287.38522541059365, 288.2, -42223.5408832171, 101320.0, 1006.0037676032673, 287.38522541059365, 288.2, -42223.5408832171, 101320.0, 1006.0037676032673, 287.38522541059365, 323.43854933808694, -6742.494176915617, 145201.140583973, 1007.9016199581774, 287.38522541059365, 323.43854933808694, -6742.494176915617, 145201.140583973, 1007.9016199581774, 287.38522541059365, 0.0, 273.900963809446, 169.59284167300237, 84791.74168172698, 1005.4614514053877, 287.38522541059365, 0.5108355114728461, 291.87110659164034, 252.14127189957713, 101320.0, 1006.1626933649725, 287.38522541059365, 0.7358330501171156, 291.87110659164034, 252.14127189957713, 101320.0, 1006.1626933649725, 287.38522541059365, 0.7358330501171156, 0.7403038111924396, 0.8903571909715106, 0.884644159087723)

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

        engine.enginecalc!(ac, "design", 1, ipcruise1, true)

        @test pare[iemfuel,ipcruise1,1] ≈ 0.00054817576594632
        @test pare[iePfan,ipcruise1,1] ≈ 6.043503411435477e6
        @test ac.engine.data.number_cells ≈ 265.5192500533146
        @test ac.engine.data.area_cell ≈ 5.0
        @test ac.engine.data.FC_heat[ipcruise1,1] ≈ 2.6919238788173264e6

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
        engine.enginecalc!(ac, "off_design", 1, iprotate, true)

        @test pare[iemfuel,iprotate,1] ≈ 0.0010447183729991173
        @test pare[iePfan,iprotate,1] ≈ 1.0000000000384096e7
        @test ac.engine.data.FC_heat[iprotate,1] ≈ 6.64805697887804e6
    end
end
