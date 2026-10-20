#----------------------------------------------------------------
#     Heat exchanger sizing tests
#
#     Thos script determines the geometric parameters for a cross-flow heat 
#     exchanger in air. It assumes that air in flowing through the exchanger
#     and heat is absorbed by a coolant fluid.

#---------------------------------     
# Test functions individually
#---------------------------------
@testset "Heat Exchanging" begin
    @testset "HEX functions" begin
        #---------------------------------     
        # hxsize!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 1144/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.Mc_in = 0.0285
        HXgas.pp_in = 40e3
        HXgas.pc_in = 1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.is_concentric = true
        HXgeom.has_recirculation = false
        HXgeom.has_shaft = false
        HXgeom.D_i = 0.564
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.n_stages = 4
        HXgeom.xt_D = 6
        HXgeom.xl_D = 1.25
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        size_out = [HXgas.Tp_out, HXgas.Tc_out, HXgas.Δp_p, HXgeom.N_t, HXgeom.n_passes, HXgeom.tD_o, HXgeom.A_cs]

        size_out_check = 
        [731.5742945899459, 666.1107129924375, 1414.5439905520996, 62.030579130833814, 8.043854538298081, 0.004760711790463513, 1.0189003226900804]

        for i = 1: length(size_out)
            @test size_out[i] ≈ size_out_check[i]
        end
        #---------------------------------     
        # hxweight()
        #---------------------------------
        gee = 9.81
        fouter = 1

        W = TASOPT.engine.hxweight(gee, HXgeom, fouter)

        W_check = 801.3484476824821

        @test W == W_check
        #---------------------------------     
        # hxoper!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 2*1144/60
        HXgas.mdot_c = 2*9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.Mc_in = 0.0285
        HXgas.pp_in = 3*40e3
        HXgas.pc_in = 3*1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.is_concentric = true
        HXgeom.has_recirculation = false
        HXgeom.D_i = 0.564
        HXgeom.t = 0.03e-2 #m, wall thicknesss
        HXgeom.tD_o = 0.004760326082769499
        HXgeom.A_cs = 1.0189779296746375
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.n_stages = 4
        HXgeom.n_passes = 8
        HXgeom.N_t = 62
        HXgeom.xt_D = 6
        HXgeom.xl_D = 1.25
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")

        TASOPT.engine.hxoper!(HXgas, HXgeom)

        oper_out = [HXgas.Tp_out, HXgas.Tc_out, HXgas.Δp_p, HXgas.ε]

        oper_out_check = [740.2990485214676, 591.3045686788635, 1706.8089965901147, 0.6502797871380385]

        for i = 1: length(oper_out)
            @test oper_out[i] ≈ oper_out_check[i]
        end

        #---------------------------------     
        # hxoptim!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 1144/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.pp_in = 40e3
        HXgas.pc_in = 1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.is_concentric = true
        HXgeom.has_recirculation = false
        HXgeom.D_i = 0.564
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.xl_D = 1
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6
        HXgeom.maxL = 0.5

        #Calculate starting point
        #First calculate minimum tube length
        _, _, _, _, cp_p_in, Rp = TASOPT.engine.gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
        γ_p_in = cp_p_in / (cp_p_in - Rp)
        ρ_p_in = HXgas.pp_in / (Rp * HXgas.Tp_in)
        Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rp * HXgas.Tp_in)

        A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream

        if HXgeom.is_concentric #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            linit = 1.1 * lmin
        else #square cross-section
            AR_min = 0.1 #Minimum aspect ratio
            lmin = sqrt(AR_min * A_cs)
            linit = sqrt(A_cs)
        end

        #Now set starting point
        initial_x = [3, 4, 4, linit] #Initial guess
        
        # Optimize heat exchanger design parameters
        TASOPT.engine.hxoptim!(HXgas, HXgeom, initial_x)

        optim_out = [HXgas.Mc_in, HXgeom.n_stages, HXgeom.xt_D, HXgeom.l]

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        Iobj = HXgas.Pl_p + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. Check I too
        I_check = 72788.03300427314
        @test Iobj ≈ I_check    rtol = 1e-5

        #---------------------------------     
        # hxoptim!() with recirculation and rectangular
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 2 * 49.9/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.825
        HXgas.Tp_in = 791
        HXgas.Tc_in = 20
        HXgas.Mp_in  = 0.01
        HXgas.pp_in = 1515e3
        HXgas.pc_in = 1515e3
        HXgas.recircT  = 200
        HXgas.h_lat  = 446e3 + 670e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.is_concentric = false
        HXgeom.has_recirculation = true
        HXgeom.t = 0.03e-2 #m, wall thicknesss
        HXgeom.xl_D = 1
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6
        HXgeom.maxL = 0.5

        #Calculate starting point
        #First calculate minimum tube length
        _, _, _, _, cp_p_in, Rh = TASOPT.engine.gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
        γ_p_in = cp_p_in / (cp_p_in - Rh)
        ρ_p_in = HXgas.pp_in / (Rh * HXgas.Tp_in)
        Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rh * HXgas.Tp_in)

        A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream

        if HXgeom.is_concentric #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            linit = 1.1 * lmin
        else #square cross-section
            AR_min = 0.1 #Minimum aspect ratio
            lmin = sqrt(AR_min * A_cs)
            linit = 1.1 * lmin
        end

        #Now set starting point
        initial_x = [2, 8, 2, linit] #Initial guess

        TASOPT.engine.hxoptim!(HXgas, HXgeom, initial_x)

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        optimrec_out = [HXgas.Mc_in, HXgas.mdot_r, HXgeom.l, HXgeom.n_stages, HXgeom.xt_D]

        Iobj_rec = HXgas.Pl_p + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. 

        I_check_rec = 3534.375943558577

        @test Iobj_rec ≈ I_check_rec    rtol = 1e-5
    end

    @testset "Thermal and pressure models" begin
        #Colburn j-factor
        Re = 1e4
        j,Cf = TASOPT.engine.jcalc_pipe(Re)

        Cf_check = 0.0791 * Re^(-0.25)
        j_check = Cf_check/2
        @test Cf ≈ Cf_check
        @test j ≈ j_check

        #Nu past tubes
        Res = [20, 500, 1e4, 3e5]
        Pr = 0.7
        N_L = 10.0
        xt_D = 3.0
        xl_D = 1.0
        Nu_check = [2.970870635447898, 13.683013469507236, 86.34877055350766, 787.9216560234382]

        for (i,Re) in enumerate(Res)
            Nu = TASOPT.engine.Nu_calc_staggered_cyl(Re, Pr, N_L, xt_D, xl_D)
            @test Nu ≈ Nu_check[i]
        end

        #Δp past tubes
        Re = 25e3
        Acs = pi * (1.265^2 - 0.564^2)/4
        l = 0.6
        tD_o = 4.78e-3
        L = 0.19
        N_tubes_tot = 1984
        xt_D = 6.0
        xl_D = 1.25
        μ_μw = 1.3

        G = 1144/60 / (Acs - 62*tD_o*l) #From Brewer
        NFV = Acs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
        Ah = N_tubes_tot * tD_o * pi * l
        Dv = 4 * NFV / Ah
        ρ = 40e3 / (287 * 750)
        
        Δp = TASOPT.engine.Δp_calc_staggered_cyl(Re, G, L, ρ, Dv, tD_o, xt_D, xl_D, μ_μw)
        Δp_check = 1323.2936826807222

        @test Δp ≈ Δp_check

        #Tube sizing
        HXgeom = TASOPT.engine.HX_tubular()
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")

        Acc = 1.372e-5 * 4 * 62
        b = pi *  0.564
        n_stages = 4
        K = pi * b * n_stages / (4 * xt_D * Acc)
        HXgeom.Δpdes = 3e6
        TASOPT.engine.tubesize!(K, HXgeom)

        t_check = 300e-6
        tDo_check = 0.004792450000885403
        @test HXgeom.t ≈ t_check
        @test HXgeom.tD_o ≈ tDo_check

        K = 10 #different case with wall thickness above minimum
        t_check = 0.0014766145463702754
        tDo_check = 0.10582404248986974

        TASOPT.engine.tubesize!(K, HXgeom)
        @test HXgeom.t ≈ t_check
        @test HXgeom.tD_o ≈ tDo_check

    end

    @testset "HEX design and off-design performance" begin
        ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
        pare = ac.pare

        pare[ieDi, :] .= 0.564
        pare[ieTft, :] .= 20
        pare[ieetab,:] .= 1.0

        ac.options.ifuel = 40
       
        # pare[iePreCMp,:] .= 0.1
        # pare[ieInterCorder,:] .= 2
        
        # pare[ieInterCMp,:] .= 0.1
        # pare[ieRegenorder,:] .= 4
        
        # pare[ieRegenMp,:] .= 0.2
        # pare[ieTurbCorder,:] .= 3
        
        # pare[ieTurbCMp,:] .= 0.02

        pare[iemcore, :] = [58.387756730737166, 59.361938832270724, 0.0, 0.0, 57.376289647792646, 48.637171035592324, 39.613606874361615, 31.882530992724078, 25.687220812619703, 22.634488608433564, 18.530429817478, 9.672011605618463, 13.622655467186979, 19.12825348500402, 24.22920161131781, 20.419257819838812, 0.0]
        pare[iemofft, :] =  [0.567, 0.567, 0.0, 0.0, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.567, 0.0]  
        pare[iefc, :] = [0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698, 0.19679538333903698] 
        pare[ieff, :] = [0.023996472039243946, 0.023927727745053126, 0.0, 0.0, 0.022840403299212705, 0.023037385319670346, 0.023344756299517953, 0.023667440439672328, 0.023928658516234116, 0.019558537480589663, 0.019290561942732494, 0.00829817423996723, 0.008666707868623658, 0.009775009482092426, 0.010029480322352872, 0.007357517728645679, 0.0]./2

        pare[ieTt19, :] = [288.2, 290.9863753875631, 0.0, 0.0, 292.16774064072325, 282.2983267937044, 268.743757672191, 256.73274363466, 247.53108063504015, 247.50937622665612, 244.44265562427674, 245.35596255379323, 255.60664560188206, 272.6173697347941, 287.3460773332474, 291.31530929825726, 0.0]
        pare[iept19, :] =  [101117.36, 104580.4483554955, 0.0, 0.0, 106074.10555835541, 84255.27417424513, 62920.64532013172, 47294.56925205369, 36374.3408691923, 36363.19484539293, 29381.955342627276, 29767.417883556787, 44295.609511469585, 64017.56898385295, 88283.93755201918, 104994.80883084564, 0.0]
        pare[ieTt25, :] = [389.95799740202244, 393.3431060406705, 0.0, 0.0, 392.78423477043805, 381.6148034607835, 364.2127150038655, 345.0356045468591, 327.66979838662934, 333.26510471570236, 328.40652447913635, 297.0142417648663, 308.4146914179435, 328.8983083692069, 341.11533726897443, 313.7329547063874, 0.0]
        pare[iept25, :] =  [257315.4699889391, 265382.5822212929, 0.0, 0.0, 264663.76651357394, 213683.54601049173, 160680.53765287113, 117608.39781182756, 86281.71218919395, 90907.98760517949, 72943.23006533818, 53595.727282542925, 78973.89139767873, 114146.15153882065, 149830.95618412812, 131942.22648489228, 0.0]
        pare[ieTt3, :] =  [850.0526625609432, 853.2327356006718, 0.0, 0.0, 840.7531226760383, 830.2224797975082, 813.7891660717617, 795.6405902172543, 779.1664744983456, 736.0498520859062, 728.3158002037491, 548.5987991211431, 566.4802195421626, 606.7099720041882, 624.817270263933, 552.6883546750172, 0.0]
        pare[iept3, :] = [3.0822488284805836e6, 3.134845596198526e6, 0.0, 0.0, 2.9881792587693473e6, 2.526639116338107e6, 2.0499960854051604e6, 1.6420310810391742e6, 1.3156760737526976e6, 1.0908958451597393e6, 884240.959052453, 355162.97269855434, 516352.83247378777, 762896.214022425, 984223.2392171419, 757513.5853096871, 0.0]
        pare[ieTt4, :] = [1833.0, 1833.0, 1833.0, 1587.0, 1783.800000500471, 1783.800000500471, 1783.800000500471, 1783.800000500471, 1783.800000500471, 1587.0000027196975, 1574.533165729207, 969.4335330427257, 995.4348851793429, 1075.9794968072274, 1100.3581190368097, 917.1708155839043, 1587.0]
        pare[ieTt49, :] =  [960.8605760851929, 961.0667107064039, 0.0, 0.0, 932.5623722599007, 931.379114632356, 930.1373020210664, 928.7863005237981, 927.5751258835901, 813.5557634539438, 806.0067604348123, 496.65626350066435, 518.3928731493768, 568.4261428950612, 599.2792291588116, 566.1504455492202, 0.0]
        pare[iept49, :] =  [187280.4402482879, 190414.81523352303, 0.0, 0.0, 181171.69203683006, 152997.44286231213, 124245.76373035634, 99618.48689089209, 79891.17968319073, 65200.92315926424, 52808.02414631647, 24694.846569804093, 37740.31793925859, 57359.1351869912, 82195.04563774749, 107666.7265230664, 0.0]
 
        ipdes = ipcruise1

        #Test precooler
        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "PreC"
        HXs[1].design_effectiveness = 0.5
        HXs[1].design_Mach = 0.1
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HX = HXs[1]

        @test HX.HXgeom.n_stages ≈ 16.13442857814795   rtol = 1e-5
        @test HX.HXgeom.n_passes ≈ 1.0725716070943603    rtol = 1e-5
        @test HX.HXgeom.l ≈ 0.4432789049197145    rtol = 1e-5
        @test HX.HXgeom.N_t ≈ 132.01600664430018    rtol = 1e-5

        @test HX.HXgas_mission[ipdes].ε ≈ 0.500000000011525    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δh_p ≈  -13469.833152449006    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δp_p ≈ 41.71181359063386   rtol = 1e-5

        for ip =1:iptotal
            @test pare[iePreCDeltah, ip] ≈ HX.HXgas_mission[ip].Δh_p
            @test pare[iePreCDeltap, ip] ≈ HX.HXgas_mission[ip].Δp_p
        end

        #Test intercooler
        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "InterC"
        HXs[1].design_effectiveness = 0.5
        HXs[1].design_Mach = 0.1
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HX = HXs[1]

        @test HX.HXgeom.n_stages ≈ 14.985782286239004     rtol = 1e-5
        @test HX.HXgeom.n_passes ≈ 1.2701493981457253    rtol = 1e-5
        @test HX.HXgeom.l ≈ 0.24878508233070346    rtol = 1e-5
        @test HX.HXgeom.N_t ≈ 148.75936468816855   rtol = 1e-5

        @test HX.HXgas_mission[ipdes].ε ≈ 0.5000000000012754    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δh_p ≈ -19023.600308918238    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δp_p ≈ 169.4607972020096     rtol = 1e-5

        for ip =1:iptotal
            @test pare[ieInterCDeltah, ip] ≈ HX.HXgas_mission[ip].Δh_p
            @test pare[ieInterCDeltap, ip] ≈ HX.HXgas_mission[ip].Δp_p
        end

        #Test cooler of turbine cool. air
        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "TurbC"
        HXs[1].design_effectiveness = 0.5
        HXs[1].design_Mach = 0.02
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HX = HXs[1]

        @test HX.HXgeom.n_stages ≈ 20.0    rtol = 1e-5
        @test HX.HXgeom.n_passes ≈  4.222920524653855    rtol = 1e-5
        @test HX.HXgeom.l ≈  0.0897318065436131    rtol = 1e-5
        @test HX.HXgeom.N_t ≈ 63.42061741899705    rtol = 1e-5

        @test HX.HXgas_mission[ipdes].ε ≈ 0.5000000000012754    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δh_p ≈ -215422.60328655195    rtol = 1e-5 
        @test HX.HXgas_mission[ipdes].Δp_p ≈ 1110.3823998901757    rtol = 1e-5

        for ip =1:iptotal
            @test pare[ieTurbCDeltah, ip] ≈ HX.HXgas_mission[ip].Δh_p
            @test pare[ieTurbCDeltap, ip] ≈ HX.HXgas_mission[ip].Δp_p
        end

        #Test regenerative cooler
        pare[ieTfuel, :] .= 20
        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "Regen"
        HXs[1].design_effectiveness = 0.5
        HXs[1].design_Mach = 0.2
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HX = HXs[1]

        @test HX.HXgeom.n_stages ≈ 6.7042669026385    rtol = 1e-5
        @test HX.HXgeom.n_passes ≈  4.1673585541263085   rtol = 1e-5 
        @test HX.HXgeom.l ≈ 0.2744402764057626    rtol = 1e-5
        @test HX.HXgeom.N_t ≈  83.75039070073112    rtol = 1e-5

        @test HX.HXgas_mission[ipdes].ε ≈ 0.5000000000012754    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δh_p ≈ -48190.134937808325    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δp_p ≈ 683.106596128081    rtol = 1e-5

        for ip =1:iptotal
            @test pare[ieRegenDeltah, ip] ≈ HX.HXgas_mission[ip].Δh_p
            @test pare[ieRegenDeltap, ip] ≈ HX.HXgas_mission[ip].Δp_p
        end

        #Test regenerative cooler with recirculation
        pare[ieTfuel, :] .= 20
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "Regen"
        HXs[1].design_effectiveness = 0.8
        HXs[1].design_Mach = 0.2
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        HXs[1].has_recirculation = true
        HXs[1].recirculation_temperature = 200.0
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        HX = HXs[1]

        @test HX.HXgeom.n_stages ≈ 19.999999999998554    rtol = 1e-5
        @test HX.HXgeom.n_passes ≈ 5.164237221491714    rtol = 1e-5
        @test HX.HXgeom.l ≈  0.27419858590558355    rtol = 1e-5
        @test HX.HXgeom.N_t ≈ 122.00395093505779   rtol = 1e-5

        @test HX.HXgas_mission[ipdes].ε ≈ 0.7999999999981817    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δh_p ≈ -88222.73967367201    rtol = 1e-5
        @test HX.HXgas_mission[ipdes].Δp_p ≈ 2642.5404384047592   rtol = 1e-5

        for ip =1:iptotal
            @test pare[ieRegenDeltah, ip] ≈ HX.HXgas_mission[ip].Δh_p
            @test pare[ieRegenDeltap, ip] ≈ HX.HXgas_mission[ip].Δp_p
        end

        #Test two HXs: intercooler and regenerative cooler
        pare[ieTfuel, :] .= 20
        
        HXs = [TASOPT.engine.make_HeatExchanger(1), TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "InterC"
        HXs[1].design_effectiveness = 0.8
        HXs[1].design_Mach = 0.1
        HXs[1].order = 1
        HXs[1].maximum_length = 0.25
        HXs[1].has_recirculation = true
        HXs[1].recirculation_temperature = 200.0

        HXs[2].type = "Regen"
        HXs[2].design_effectiveness = 0.8
        HXs[2].design_Mach = 0.2
        HXs[2].order = 2
        HXs[2].maximum_length = 0.25
        ac.engine.heat_exchangers = HXs

        HXs = TASOPT.hxdesign!(ac, ipdes, 1)

        @test pare[ieInterCDeltah,:] ≈ [-47529.983813768704, -47872.611440618784, 0.0, 0.0, -46352.290847166885, -45737.52211589867, -44368.27724079696, -42610.805067712456, -40892.32336811514, -35751.99236363902, -35106.687590422494, -14678.390016351703, -15721.335516107116, -18787.94106567447, -19945.59089965135, -13489.800142848635, 0.0]
        @test pare[ieInterCDeltap,:] ≈  [6956.456335392031, 7024.684262537262, 0.0, 0.0, 6615.539945295114, 5816.831645336263, 4991.567831102287, 4263.954011671673, 3646.2867476659444, 2804.477506629798, 2348.9135713448195, 841.2918318651618, 1159.9613172683567, 1647.7867641224643, 2048.338254893675, 1547.0833751982325, 0.0]
        @test pare[ieRegenDeltah,:] ≈ [-56867.32798647176, -56102.26491999626, 0.0, 0.0, -53153.37916611304, -58273.96779443405, -65296.91476474004, -73042.86637903785, -80589.45461977209, -58633.32913037343, -60147.96557543095, -12759.553185814904, -13754.009796559461, -17240.721435134066, -18617.047029681795, -13702.598399612543, 0.0]
        @test pare[ieRegenDeltap,:] ≈ [5426.266142447696, 5507.868187373721, 0.0, 0.0, 5245.4980403172995, 4538.611990755371, 3782.0265080148633, 3115.8666195008805, 2567.8077651877375, 2136.766563284311, 1779.1545982077182, 634.5748157074588, 849.8901765252981, 1189.389772301055, 1380.5925230574924, 712.9106725171165, 0.0]

    end

    @testset "Radiator design and off-design performance" begin

        ac = load_default_model()
        pare = ac.pare
        pare[ieTt21,:,1] .= [343.5937504775422, 345.1223581889431, 0.0, 0.0, 313.78565028844844, 307.1133732298289, 299.4976209320901, 295.77161933759936, 297.7979970527857, 292.4179898692317, 290.77602320035004, 263.0449635812037, 281.20843068925603, 298.112881341098, 309.3973072464574, 298.7327491186001, 0.0]
        pare[iept21,:,1] .= [175101.1544126459, 178399.39021402935, 0.0, 0.0, 133748.37095494103, 110626.50581060312, 89022.3186448056, 73988.05771998747, 64580.336522473954, 61271.983314487086, 60150.59571341334, 44236.17815241701, 63955.92953145897, 88974.3688113112, 114550.16697333875, 114113.91098410703, 0.0]
        pare[ieRadiatorCoolantT,:,1] .= [453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15, 453.15]
        pare[ieRadiatorCoolantP,:,1] .= [300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0, 300000.0]
        pare[iemfan,:,1] .= [358.3434869671312, 366.38016489316016, 0.0, 0.0, 249.19699475758492, 234.48400704476614, 200.65335506151317, 170.25878321648813, 146.90816150458102, 141.16684807770758, 139.1072989341617, 108.25761644603271, 151.1038255310103, 202.63153614761245, 247.82443427439475, 149.89424326937296, 0.0]
        pare[ieRadiatorHeat,:,1] .= [1.3296113956684684e7, 1.3296113956683043e7, 0.0, 0.0, 1.6765644717054842e6, 1.863898537598179e6, 2.0218038992171476e6, 2.239585517913293e6, 2.592509216389226e6, 2.0910313931780618e6, 1.954040175585774e6, 379399.2116786214, 657833.3780815811, 1.0405791722436543e6, 1.2058361829043678e6, 200019.69516559975, 0.0]
        
        pare[ieDi,:,1] .= 0.4 #Inner diameter of HEX 

        HXs = [TASOPT.engine.make_HeatExchanger(1)]
        HXs[1].type = "Radiator"
        HXs[1].design_effectiveness = 0.7
        HXs[1].design_Mach = 0.12
        HXs[1].order = 1
        HXs[1].maximum_length = 2.0
        ac.engine.heat_exchangers = HXs
        HXs = TASOPT.hxdesign!(ac, ipstatic, 1)

        RadDeltah = [37104.38291823695, 36290.485213795124, 0.0, 0.0, 6727.86793972544, 7948.936735981042, 10076.103131180309, 13154.008713110605, 17647.142199845595, 14812.481978963075, 14046.99962228872, 3504.596019511984, 4353.519017601429, 5135.326869780114, 4865.687221015694, 1334.405450155591, 0.0]
        RadDeltap = [10323.181606429736, 10592.486241518414, 0.0, 0.0, 5937.103176636959, 6281.141626491474, 5716.324620345804, 5036.269927471252, 4460.559580576974, 4258.477622450499, 4189.339004377093, 3127.789243974919, 4331.454425196778, 5730.505091121955, 6729.144232579738, 2533.1020698654006, 0.0]

        for ip =1:iptotal
            @test pare[ieRadiatorDeltah, ip, 1] ≈ RadDeltah[ip]
            @test pare[ieRadiatorDeltap, ip, 1] ≈ RadDeltap[ip]
        end

    end

end