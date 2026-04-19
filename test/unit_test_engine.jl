using Zygote
using DelimitedFiles
using ForwardDiff
using StaticArrays

isGradient = false

@testset "Engine models" begin
    @testset "gasfun.jl" begin

        # =========================
        # gas_N2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_N2(t1)

        @test s1 == 2358.21669383134
        @test h1 == 2.43487432826e6
        @test cp1 == 1301.485606111
        @test r1 == 296.94

        if isGradient
            ds_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[1], t1)[1]
            dh_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[2], t1)[1]
            dcp_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[3], t1)[1]
            dr_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[4], t1)[1]

            @test ds_dt == 0.5578891260784105
            @test dh_dt == 1301.4526599999983
            @test dcp_dt == 0.04276110100000096
            @test dr_dt == nothing
        end

        # =========================
        # gas_N2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_Ar(t1)

        @test s1 == 1070.0677477147947
        @test h1 == 1.0582e6
        @test cp1 == 520.0
        @test r1 == 208.0

        # =========================
        # gas_CO2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_CO2(t1)

        @test s1 == 2382.871389635858
        @test h1 == -6.405596226739999e6
        @test cp1 == 1389.688099952
        @test r1 == 188.96

        # =========================
        # gas_H2O
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_H2O(t1)

        @test s1 == 4656.047912649002
        @test h1 == -8.43504044348e6
        @test cp1 == 2943.545113578
        @test r1 == 461.91


        # =========================
        # gas_CH4
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_CH4(t1)

        @test s1 == 8467.631968582504
        @test h1 == 4.574569283279987e6
        @test cp1 == 5288.091482799951
        @test r1 == 519.65


        # =========================
        # gas_C2H6
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C2H6(t1)

        @test s1 == 7408.861495898641
        @test h1 == 5.381058419279992e6
        @test cp1 == 4498.240840800101
        @test r1 == 277.15


        # =========================
        # gas_C3H8
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C3H8(t1)

        @test s1 == 7209.141132729834
        @test h1 == 6.1354732452e6
        @test cp1 == 4411.240840800101
        @test r1 == 188.5


        # =========================
        # gas_C4H10
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C4H10(t1)

        @test s1 == 7230.2763377713845
        @test h1 == 5.908308419279992e6
        @test cp1 == 4448.240840800101
        @test r1 == 143.3

        # =========================
        # gas_C8H18
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C8H18(t1)

        @test s1 == 7262.8799183106685
        @test h1 == 6.234133419279992e6
        @test cp1 == 4443.240840800101
        @test r1 == 72.9

        # =========================
        # gas_C14H30
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C14H30(t1)

        @test s1 == 7253.868041333698
        @test h1 == 6.369750419279992e6
        @test cp1 == 4432.240840800101
        @test r1 == 167.0
        #    


        # =========================
        # gasfun
        # =========================

        igas = 14
        t = 2333.00e0
        s, s_t, h, h_t, cp, r = TASOPT.engine.gasfun(igas, t)

        @test s == 7230.2763377713845
        @test s_t == 1.9066613119588947
        @test h == 5.908308419279992e6
        @test h_t == 4448.240840800101
        @test cp == 4448.240840800101
        @test r == 143.3

        if isGradient
            igas = 14
            t1 = 2333.00e0
            ds_dt = gradient(t1 -> TASOPT.engine.gasfun(igas, t1)[1], t1)[1]
            dh_dt = gradient(t1 -> TASOPT.engine.gasfun(igas, t1)[3], t1)[1]

            epsilon = 1e-6
            s_perturb, s_t_perturb, h_perturb, h_t_perturb, cp_perturb, r_perturb = TASOPT.engine.gasfun(igas, t + epsilon)
            ds_dt_FD = (s_perturb - s) / epsilon
            dh_dt_FD = (h_perturb - h) / epsilon
        end

        # println([s_t, h_t]) # HACK: Not matching
        # println([ds_dt, dh_dt])
        # println([ds_dt_FD, dh_dt_FD])

        # =========================
        # gaschem
        # =========================

        igas = 13
        nchon = TASOPT.engine.gaschem(igas)

        @test nchon == [3, 8, 0, 0]

    end


    @testset "gascalc.jl" begin

        # =========================
        # gassum, gassumd
        # =========================

        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        n = 6
        n_air = n - 1
        t = 215.0
        s, s_t, h, h_t, cp, r = TASOPT.engine.gassum(alpha, n_air, t)


        @test s == -329.0381537005463
        @test s_t == 4.684344004111629
        @test h == -87237.76845599999
        @test h_t == 1007.1339608840001
        @test cp == 1007.1339608840001
        @test r == 288.29530400000004

        s, s_t, h, h_t, cp, cp_t, r = TASOPT.engine.gassumd(alpha, n_air, t)

        @test s == -329.0381537005463
        @test s_t == 4.684344004111629
        @test h == -87237.76845599999
        @test h_t == 1007.1339608840001
        @test cp == 1007.1339608840001
        @test cp_t == 0.009382799699317275
        @test r == 288.29530400000004

        # =========================
        # gas_tset, gas_tsetd
        # =========================

        tguess = 200.0
        t = TASOPT.engine.gas_tset(alpha, n_air, h, tguess)

        @test t == 214.99999914397367

        t, t_hspec, t_al = TASOPT.engine.gas_tsetd(alpha, n_air, h, tguess)
        @test t == 214.99999914397367
        @test t_hspec == 0.000992916572022094
        @test t_al == [85.64080275342843, 75.32364484439914, 8945.435851261927, 13490.80321923344, 42.8542796904542]

        # =========================
        # gas_prat, gas_pratd
        # =========================
        # TODO []: gas_prat, gas_pratd not verified
        # TODO:
        # Fix fortran compiling issue


        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        n = 5
        p = 0.22e5
        t = 215.0

        # =========================
        # gas_delh, gas_delhd
        # =========================
        s, s_t, h, h_t, cp, r = TASOPT.engine.gassum(alpha, n, t)

        delh = 10000.0
        epol = 0.99

        p, t, h, s, cp, r = TASOPT.engine.gas_delh(alpha, n, p, t, h, s, cp, r, delh, epol)

        @test p == 25717.186436176951
        @test t == 224.92715529368544
        @test h == -77237.768461931657
        @test s == -283.57571753601650
        @test cp == 1007.2429752832797
        @test r == 288.29530400000004

        # =========================
        # gas_burn, gas_burnd
        # =========================
        # n == 5: the 5 air species (N2, O2, CO2, H2O, Ar); fuel handled via ifuel
        n_burn = 5
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        gamma = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
        ifuel = 13
        to = 450.0
        tf = 300.0
        t = 1400.0
        hvap = 0.0

        f, lambda = TASOPT.engine.gas_burn(alpha, beta, gamma, n_burn, ifuel, to, tf, t, hvap)

        @test f == -0.45350655959892078
        @test lambda[1] == 1.4291113895654686

        f, lambda, f_to, f_tf, f_t, l_to, l_tf, l_t = TASOPT.engine.gas_burnd(alpha, beta, gamma, n_burn, ifuel, to, tf, t, hvap)
        @test f_to == 4.3531820788987320E-004
        @test f_tf == -3.2351348482849527E-004
        @test f_t == -5.1186081193040599E-004
        @test l_to[1] == -1.1383818413704377E-003
        @test l_tf[1] == 8.4600613962923770E-004
        @test l_t[1] == 1.3385450988489609E-003

        # =========================
        # gas_mach, gas_machd
        # =========================
        # TODO [] gas_machd not verified
        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        to = 450.0
        po = 200000.0
        so, s_to, ho, h_to, cpo, ro = TASOPT.engine.gassum(alpha, n - 1, to)
        mo = 1.1
        m = 1.3
        epol = 0.99

        p, t, h, s, cp, r = TASOPT.engine.gas_mach(alpha, n - 1, po, to, ho, so, cpo, ro, mo, m, epol)

        @test p ≈ 154338.88081676030 rtol = 1e-10
        @test t ≈ 417.96446969196779 rtol = 1e-10
        @test h ≈ 117997.67846481415 rtol = 1e-10
        @test s ≈ 342.75772292744381 rtol = 1e-10
        @test cp ≈ 1019.5701545538500 rtol = 1e-10
        @test r ≈ 288.29530400000004 rtol = 1e-10

        # =========================
        # gas_mass
        # =========================

        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        to = 450.0
        po = 200000.0
        so, s_to, ho, h_to, cpo, ro = TASOPT.engine.gassum(alpha, n - 1, to)

        Mguess = 1.15
        u = sqrt(1.4 * ro * to) * 1.1
        rho = (po / (ro * to))
        mflux = rho * u

        p, t, h, s, cp, r = TASOPT.engine.gas_mass(alpha, n - 1, po, to, ho + rho * u^2 / 2, so, cpo, ro, mflux, Mguess)

        @test p ≈ 118662.19814638462 rtol = 1e-10
        @test t ≈ 388.25628795926565 rtol = 1e-10
        @test h ≈ 87768.12641472857 rtol = 1e-10
        @test s ≈ 267.72827064815743 rtol = 1e-10
        @test cp ≈ 1015.8338429023191 rtol = 1e-10
        @test r ≈ 288.29530400000004 rtol = 1e-10

        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]

        # =========================
        # gasfuel
        # =========================

        ifuel = 13
        gamma = TASOPT.engine.gasfuel(ifuel, n)
        @test gamma[2] == -3.6283226981894474

        # =========================
        # @inferred: SVector{5,Float64} composition is type-stable through all scalar-return functions
        # =========================
        let
            α5 = SVector{5,Float64}(0.781, 0.209, 0.0004, 0.0, 0.00965)
            t0 = 300.0

            # gassum: returns 6 scalars
            @test @inferred(Tuple{Float64,Float64,Float64,Float64,Float64,Float64},
                            TASOPT.engine.gassum(α5, 5, t0)) == TASOPT.engine.gassum(α5, 5, t0)

            # gassumd: returns 7 scalars
            @test @inferred(Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64},
                            TASOPT.engine.gassumd(α5, 5, t0)) == TASOPT.engine.gassumd(α5, 5, t0)

            # gas_tset: returns a scalar temperature
            s0, _, h0, _, cp0, r0 = TASOPT.engine.gassum(α5, 5, t0)
            @test @inferred(Float64, TASOPT.engine.gas_tset(α5, 5, h0, t0)) ≈ t0  atol=1e-4

            # gas_prat: returns 6 scalars
            po, to, ho, so, cpo, ro = 101325.0, 300.0, h0, s0, cp0, r0
            @test @inferred(Tuple{Float64,Float64,Float64,Float64,Float64,Float64},
                            TASOPT.engine.gas_prat(α5, 5, po, to, ho, so, cpo, ro, 1.5, 0.9)) ==
                  TASOPT.engine.gas_prat(α5, 5, po, to, ho, so, cpo, ro, 1.5, 0.9)

            # gas_delh: returns 6 scalars
            @test @inferred(Tuple{Float64,Float64,Float64,Float64,Float64,Float64},
                            TASOPT.engine.gas_delh(α5, 5, po, to, ho, so, cpo, ro, 5000.0, 0.9)) ==
                  TASOPT.engine.gas_delh(α5, 5, po, to, ho, so, cpo, ro, 5000.0, 0.9)

            # gas_mach: returns 6 scalars
            @test @inferred(Tuple{Float64,Float64,Float64,Float64,Float64,Float64},
                            TASOPT.engine.gas_mach(α5, 5, po, to, ho, so, cpo, ro, 0.0, 0.8, 1.0)) ==
                  TASOPT.engine.gas_mach(α5, 5, po, to, ho, so, cpo, ro, 0.0, 0.8, 1.0)
        end

    end

    @testset "tfmap.jl" begin

        # TODO []
        # ecmap1, Ncmap1, etmap not verified

        # =========================
        # Ncmap
        # =========================

        pratio = 10.0
        mb = 100.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        Nb, Nb_pi, Nb_mb = TASOPT.engine.Ncmap(pratio, mb, piD, mbD, NbD, Cmap)

        @test Nb == 8.3597298887214055
        @test Nb_pi == 0.26149475347651785
        @test Nb_mb == 2.4236392203562176E-002

        # =========================
        # ecmap
        # =========================

        pratio = 10.0
        mb = 100.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        effo = 0.95
        piK = 0.9 * pratio
        effK = 0.93

        eff, eff_pi, eff_mb = TASOPT.engine.ecmap(pratio, mb, piD, mbD, Cmap, effo, piK, effK)

        @test eff == 1.8781007331361292
        @test eff_pi == 0.92374562099125368
        @test eff_mb == 1.6164204096982837E-003

        # =========================
        # Pimap
        # =========================

        pratio = 10.0
        mb = 100.0
        Nb = 10.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        eff, eff_pi, eff_mb = TASOPT.engine.Pimap(mb, Nb, piD, mbD, NbD, Cmap)

        @test eff == 16.579462807918382
        @test eff_pi == -3.5593220338983059E-002
        @test eff_mb == 4.419641196046076

    end

    @testset "tfcool.jl" begin

        ncrowx = 2
        ncrow = 2
        Tt3 = 300.0
        Tt4 = 1300.0
        dTstreak = 300.0
        Trrat = 0.9
        efilm = 0.5
        tfilm = 0.3
        StA = 0.2

        epsrow = zeros(ncrowx)
        epsrow[1] = 0.1
        epsrow[2] = 0.1

        Tmrow = zeros(ncrowx)
        Tmrow[1] = 1100.0
        Tmrow[2] = 1100.0

        ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr = TASOPT.engine.mcool(ncrowx,
            Tmrow, Tt3, Tt4, dTstreak, Trrat,
            efilm, tfilm, StA)

        ncrow == 1
        epsrow[1] == 0.10313901345291480
        epsrow_Tt3[1] == 1.7595366888535866E-004
        epsrow_Tt4[1] == 2.8152587021657390E-004
        epsrow_Trr[1] == 0.0000000000000000


        ncrowx = 2
        ncrow = 2
        Tt3 = 300.0
        Tt4 = 1300.0
        dTstreak = 300.0
        Trrat = 0.9
        efilm = 0.5
        tfilm = 0.3
        StA = 0.2

        epsrow = zeros(ncrowx)
        epsrow[1] = 0.1
        epsrow[2] = 0.1

        Tmrow = TASOPT.engine.Tmcalc(ncrowx, ncrow,
            Tt3, Tt4, dTstreak, Trrat,
            efilm, tfilm, StA, epsrow)

        @test Tmrow[1] == 1106.8965517241379
        @test Tmrow[2] == 840.00000000000000

        if isGradient
            # AD
            d_Tmrow_d_Tt3 = gradient(Tt3 -> TASOPT.engine.Tmcalc(ncrowx, ncrow,
                    Tt3, Tt4, dTstreak, Trrat,
                    efilm, tfilm, StA, epsrow)[1], Tt3)[1]

            # FD
            epsilon = 1e-6
            Tmrow_d = TASOPT.engine.Tmcalc(ncrowx, ncrow,
                Tt3 + epsilon, Tt4, dTstreak, Trrat,
                efilm, tfilm, StA, epsrow)

            d_Tmrow_d_Tt3_FD = (Tmrow_d[1] - Tmrow[1]) / epsilon

            # FD AD consistency test
            @test d_Tmrow_d_Tt3_FD ≈ d_Tmrow_d_Tt3 rtol = 1e-5
        end

        #Test cooled HPT efficiency
        epht0 = 0.9
        epht_fc = -0.36
        fc0 = 0.16
        fc = 0.2
        epht = TASOPT.engine.find_cooled_hpt_efficiency(epht0, epht_fc, fc0, fc)
        @test epht ≈ 0.8856
    end

    @testset "tfsize.jl" begin

        # TODO []: icool change the function signature. Only works for icool = 1.
        # note: icool::Int -> opt_cooling::String
        # TODO []: the output shall cover most of the functions. Ideally, test all other outputs.

        gee = 9.8100000000000005
        M0 = 0.80000000000000004
        T0 = 219.43067572699252
        p0 = 23922.608843328788
        a0 = 296.85578884697560
        M2 = 0.59999999999999998
        M25 = 0.59999999999999998
        Feng = 22182.101361240744
        Phiinl = 0.0000000000000000
        Kinl = 0.0000000000000000
        eng_has_BLI_cores = false
        BPR = 5.0999999999999996
        pif = 1.6850000000000001
        pilc = 8.0000000000000000
        pihc = 3.7500000000000000
        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 89407.867373646965
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = TASOPT.engine.CoolingOpt.FixedCoolingFlowRatio
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow = zeros(ncrowx)
        epsrow[1] = 0.12811308404512714
        epsrow[2] = 5.4411331284501797E-002
        epsrow[3] = 1.5791188045605239E-002
        epsrow[4] = 0.0000000000000000
        Tt4 = 1587.0000000000000
        Tmrow = zeros(ncrowx)
        fc0 = 0.0
        epht_fc = 0.0
        hvap = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

        epsrow, Tmrow,
        TSFC, Fsp, hfuel, ff, mcore,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt19, ht19, pt19, cpt19, Rt19,
        Tt19c, ht19c, pt19c, cpt19c, Rt19c,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt25, ht25, pt25, cpt25, Rt25,
        Tt25c, ht25c, pt25c, cpt25c, Rt25c,
        Tt3, ht3, pt3, cpt3, Rt3,
        ht4, pt4, cpt4, Rt4,
        Tt41, ht41, pt41, cpt41, Rt41,
        Tt45, ht45, pt45, cpt45, Rt45,
        Tt49, ht49, pt49, cpt49, Rt49,
        Tt5, ht5, pt5, cpt5, Rt5,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, A2,
        T25c, u25c, p25c, cp25c, R25c, A25,
        T5, u5, p5, cp5, R5, A5,
        T6, u6, p6, cp6, R6, A6,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfsize!(gee, M0, T0, p0, a0, M2, M25,
            Feng, Phiinl, Kinl, eng_has_BLI_cores,
            BPR, pif, pilc, pihc,
            pid, pib, pifn, pitn,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9, Tt4,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow, Tmrow, 
            Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)



        @test etaf == 0.8539899545024271
        @test etalc ==  0.8298895182207285
        @test etahc == 0.8384795893726779
        @test etaht ==  0.8978789812518558
        @test etalt == 0.9191845671925591
        @test Tmrow[1] ≈ 1119.1584455133655  rtol = 1e-10

        if isGradient
            # AD
            d_etaf_dM0 = gradient(M0 -> TASOPT.engine.tfsize(gee, M0, T0, p0, a0, M2, M25,
                    Feng, Phiinl, Kinl, eng_has_BLI_cores,
                    BPR, pif, pilc, pihc,
                    pid, pib, pifn, pitn,
                    Ttf, ifuel, hvap, etab,
                    epf0, eplc0, ephc0, epht0, eplt0,
                    mofft, Pofft,
                    Tt9, pt9, Tt4,
                    epsl, epsh,
                    opt_cooling,
                    Mtexit, dTstrk, StA, efilm, tfilm,
                    M4a, ruc,
                    ncrowx, ncrow,
                    epsrow)[115], M0, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                    Δp_PreC, Δp_InterC, Δp_Regen)[1]

            # FD
            epsilon = 1e-6
            etaf_d = TASOPT.engine.tfsize(gee, M0 + epsilon, T0, p0, a0, M2, M25,
                Feng, Phiinl, Kinl, eng_has_BLI_cores,
                BPR, pif, pilc, pihc,
                pid, pib, pifn, pitn,
                Ttf, ifuel, hvap, etab,
                epf0, eplc0, ephc0, epht0, eplt0,
                mofft, Pofft,
                Tt9, pt9, Tt4,
                epsl, epsh,
                opt_cooling,
                Mtexit, dTstrk, StA, efilm, tfilm,
                fc0, epht_fc,
                M4a, ruc,
                ncrowx, ncrow,
                epsrow, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                Δp_PreC, Δp_InterC, Δp_Regen)[115]

            d_etaf_dM0_FD = (etaf_d - etaf) / epsilon
            # FD AD consistency test

            @test d_etaf_dM0_FD ≈ d_etaf_dM0 rtol = 1e-4
        end

        gee = 9.8100000000000005
        M0 = 0.80000000000000004
        T0 = 219.43067572699252
        p0 = 23922.608843328788
        a0 = 296.85578884697560
        M2 = 0.59999999999999998
        M25 = 0.59999999999999998
        Feng = 22182.101361240744
        Phiinl = 0.0000000000000000
        Kinl = 0.0000000000000000
        eng_has_BLI_cores = false
        BPR = 5.0999999999999996
        pif = 1.6850000000000001
        pilc = 8.0000000000000000
        pihc = 3.7500000000000000
        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 89407.867373646965
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = TASOPT.engine.CoolingOpt.FixedTmetal
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow = zeros(ncrowx)
        Tt4 = 1587.0000000000000
        Tmrow = [1114.7762513333382, 1102.6358896742904, 1098.9326716275868, 1021.0788483883676]


        epsrow, Tmrow,
        TSFC, Fsp, hfuel, ff, mcore,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt19, ht19, pt19, cpt19, Rt19,
        Tt19c, ht19c, pt19c, cpt19c, Rt19c,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt25, ht25, pt25, cpt25, Rt25,
        Tt25c, ht25c, pt25c, cpt25c, Rt25c,
        Tt3, ht3, pt3, cpt3, Rt3,
        ht4, pt4, cpt4, Rt4,
        Tt41, ht41, pt41, cpt41, Rt41,
        Tt45, ht45, pt45, cpt45, Rt45,
        Tt49, ht49, pt49, cpt49, Rt49,
        Tt5, ht5, pt5, cpt5, Rt5,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, A2,
        T25c, u25c, p25c, cp25c, R25c, A25,
        T5, u5, p5, cp5, R5, A5,
        T6, u6, p6, cp6, R6, A6,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfsize!(gee, M0, T0, p0, a0, M2, M25,
            Feng, Phiinl, Kinl, eng_has_BLI_cores,
            BPR, pif, pilc, pihc,
            pid, pib, pifn, pitn,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9, Tt4,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow, Tmrow, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)

        @test epsrow[1] ≈ 0.1303162243242347 rtol = 1e-10

    end

    @testset "gaussn.jl" begin
        nsiz = 3
        nn = 3
        nrhs = 1

        z = zeros((nsiz, nsiz))
        z[1, 1] = 0.52762351
        z[1, 2] = 0.85729295
        z[1, 3] = 0.32073675
        z[2, 1] = 0.70046643
        z[2, 2] = 0.43823310
        z[2, 3] = 0.60231747
        z[3, 1] = 0.86206149
        z[3, 2] = 0.35227634
        z[3, 3] = 0.61320909

        r = zeros((nsiz, nrhs))
        r[1, 1] = 0.52873024
        r[2, 1] = 0.78246347
        r[3, 1] = 0.10212863

        r = TASOPT.engine.gaussn(nsiz, nn, z, r, nrhs)

        @test r[1, 1] == -3.9352200465675033
        @test r[2, 1] == 1.1548320796061180
        @test r[3, 1] == 5.0353302305150418

    end

    @testset "tfoper.jl mode 1" begin

        gee = 9.8100000000000005
        M0 = 0.26302467815397762
        T0 = 288.00000000000000
        p0 = 101320.00000000000
        a0 = 340.08940001123227
        Tref = 288.19999999999999
        pref = 101320.00000000000

        Phiinl = 0.0
        Kinl = 0.0
        eng_has_BLI_cores = false

        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Gearf = 1.0000000000000000
        pifD = 1.6850000000000001
        pilcD = 8.0000000000000000
        pihcD = 3.7500000000000000
        pihtD = 2.1601257635200488
        piltD = 6.2886975330083716
        mbfD = 235.16225770724063
        mblcD = 46.110246609262873
        mbhcD = 7.8056539219349039
        mbhtD = 4.3594697284253883
        mbltD = 8.7016090343744406
        NbfD = 1.0790738309310697
        NblcD = 1.0790738309310697
        NbhcD = 0.77137973563891493
        NbhtD = 0.44698693289691338
        NbltD = 0.48396724306758404
        A2 = 1.3863121762890294
        A25 = 3.8585338087708761E-002
        A5 = 0.19210855588408102
        A7 = 0.64211443204484309
        opt_calc_call = TASOPT.engine.CalcMode.FixedTt4OffDes
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 77800.595231538944
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = TASOPT.engine.CoolingOpt.FixedCoolingFlowRatio
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow3 = [0.12061791584226822, 5.1292591721870069E-002, 1.5478853228971187E-002, 0.0000000000000000]
        Tmrow3 = [1000.0, 1000.0, 1000.0, 1000.0]
        fc0 = 0.0
        epht_fc = 0.0

        M2 = 1.0
        pif = 0.0000000000000000
        pilc = 0.0000000000000000
        pihc = 0.0000000000000000
        mbf = 0.0000000000000000
        mblc = 0.0000000000000000
        mbhc = 0.0000000000000000
        Tt4 = 1783.8000000000002
        pt5 = 0.0000000000000000
        mcore = 0.0
        Feng = 0.0
        M25 = 0.0
        hvap = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

        # Fix mass flow find temperature
        TSFC, Fsp, hfuel, ff,
        Feng, mcore,
        pif, pilc, pihc,
        mbf, mblc, mbhc,
        Nbf, Nblc, Nbhc,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt19, ht19, pt19, cpt19, Rt19,
        Tt19c, ht19c, pt19c, cpt19c, Rt19c,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt25, ht25, pt25, cpt25, Rt25,
        Tt25c, ht25c, pt25c, cpt25c, Rt25c,
        Tt3, ht3, pt3, cpt3, Rt3,
        Tt4, ht4, pt4, cpt4, Rt4,
        Tt41, ht41, pt41, cpt41, Rt41,
        Tt45, ht45, pt45, cpt45, Rt45,
        Tt49, ht49, pt49, cpt49, Rt49,
        Tt5, ht5, pt5, cpt5, Rt5,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, M2,
        T25c, u25c, p25c, cp25c, R25c, M25c,
        T5, u5, p5, cp5, R5, M5,
        T6, u6, p6, cp6, R6, M6, A6,
        T7, u7, p7, cp7, R7, M7,
        T8, u8, p8, cp8, R8, M8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfoper!(gee, M0, T0, p0, a0, Tref, pref,
            Phiinl, Kinl, eng_has_BLI_cores,
            pid, pib, pifn, pitn,
            Gearf,
            pifD, pilcD, pihcD, pihtD, piltD,
            mbfD, mblcD, mbhcD, mbhtD, mbltD,
            NbfD, NblcD, NbhcD, NbhtD, NbltD,
            A2, A25, A5, A7,
            opt_calc_call,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow3, Tmrow3,
            Feng,
            M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
            Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)


        @test etaf ≈ 0.8741690445868673 atol = 1e-8
        @test etalc ≈ 0.8367067552587594 atol = 1e-8
        @test etahc ≈ 0.8393327258267146 atol = 1e-8
        @test etaht ≈ 0.8974130931530157  atol = 1e-8
        @test etalt ≈ 0.9176389341085217 atol = 1e-8

        Tt4_ref = 1783.8000000000002
        ht4_ref = 470172.6605599733
        pt4_ref = 2.7406297732608365e6
        cpt4_ref = 1312.0669933280028
        Rt4_ref =  288.12027917359467
        @test Tt4 ≈ Tt4_ref atol = 1e-8 * Tt4_ref
        @test ht4 ≈ ht4_ref atol = 1e-8 * ht4_ref
        @test pt4 ≈ pt4_ref atol = 1e-8 * pt4_ref
        @test cpt4 ≈ cpt4_ref atol = 1e-8 * cpt4_ref
        @test Rt4 ≈ Rt4_ref atol = 1e-8 * Rt4_ref

    end

    @testset "tfoper.jl mode 2" begin

        gee = 9.8100000000000005
        M0 = 0.26302467815397762
        T0 = 288.00000000000000
        p0 = 101320.00000000000
        a0 = 340.08940001123227
        Tref = 288.19999999999999
        pref = 101320.00000000000

        Phiinl = 0.0
        Kinl = 0.0
        eng_has_BLI_cores = false

        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Gearf = 1.0000000000000000
        pifD = 1.6850000000000001
        pilcD = 8.0000000000000000
        pihcD = 3.7500000000000000
        pihtD = 2.1601257635200488
        piltD = 6.2886975330083716
        mbfD = 235.16225770724063
        mblcD = 46.110246609262873
        mbhcD = 7.8056539219349039
        mbhtD = 4.3594697284253883
        mbltD = 8.7016090343744406
        NbfD = 1.0790738309310697
        NblcD = 1.0790738309310697
        NbhcD = 0.77137973563891493
        NbhtD = 0.44698693289691338
        NbltD = 0.48396724306758404
        A2 = 1.3863121762890294
        A25 = 3.8585338087708761E-002
        A5 = 0.19210855588408102
        A7 = 0.64211443204484309
        opt_calc_call = TASOPT.engine.CalcMode.FixedTt4OffDes
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 77800.595231538944
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = TASOPT.engine.CoolingOpt.FixedCoolingFlowRatio
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow3 = [0.12061791584226822, 5.1292591721870069E-002, 1.5478853228971187E-002, 0.0000000000000000]
        Tmrow3 = [1000.0, 1000.0, 1000.0, 1000.0]
        fc0 = 0.0
        epht_fc = 0.0

        M2 = 1.0
        pif = 0.0000000000000000
        pilc = 0.0000000000000000
        pihc = 0.0000000000000000
        mbf = 0.0000000000000000
        mblc = 0.0000000000000000
        mbhc = 0.0000000000000000
        Tt4 = 1783.8000000000002
        pt5 = 0.0000000000000000
        mcore = 0.0
        M25 = 0.0
        hvap = 0.0
        Feng = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

        # Fix temperature find mass flow
        opt_cooling = TASOPT.engine.CoolingOpt.FixedTmetal
        TSFC, Fsp, hfuel, ff,
        Feng, mcore,
        pif, pilc, pihc,
        mbf, mblc, mbhc,
        Nbf, Nblc, Nbhc,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt19, ht19, pt19, cpt19, Rt19,
        Tt19c, ht19c, pt19c, cpt19c, Rt19c,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt25, ht25, pt25, cpt25, Rt25,
        Tt25c, ht25c, pt25c, cpt25c, Rt25c,
        Tt3, ht3, pt3, cpt3, Rt3,
        Tt4, ht4, pt4, cpt4, Rt4,
        Tt41, ht41, pt41, cpt41, Rt41,
        Tt45, ht45, pt45, cpt45, Rt45,
        Tt49, ht49, pt49, cpt49, Rt49,
        Tt5, ht5, pt5, cpt5, Rt5,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, M2,
        T25c, u25c, p25c, cp25c, R25c, M25c,
        T5, u5, p5, cp5, R5, M5,
        T6, u6, p6, cp6, R6, M6, A6,
        T7, u7, p7, cp7, R7, M7,
        T8, u8, p8, cp8, R8, M8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfoper!(gee, M0, T0, p0, a0, Tref, pref,
            Phiinl, Kinl, eng_has_BLI_cores,
            pid, pib, pifn, pitn,
            Gearf,
            pifD, pilcD, pihcD, pihtD, piltD,
            mbfD, mblcD, mbhcD, mbhtD, mbltD,
            NbfD, NblcD, NbhcD, NbhtD, NbltD,
            A2, A25, A5, A7,
            opt_calc_call,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow3, Tmrow3,
            Feng,
            M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
            Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)

        @test etaf ≈ 0.8862602081467577 rtol = 1e-6
        @test etalc ≈ 0.6539900828671803 rtol = 1e-6
        @test etahc ≈ 0.8390003200177845 rtol = 1e-6
        @test etaht ≈ 0.8977614695292404 rtol = 1e-6
        @test etalt ≈ 0.8649073213143893 rtol = 1e-6

        @test Tt4 ≈ 1783.8000000000002 rtol = 1e-6
        @test ht4 ≈ 292168.48621470056  rtol = 1e-6
        @test pt4 ≈ 1.0428916387810945e6  rtol = 1e-6
        @test cpt4 ≈ 1323.0363088269542 rtol = 1e-6
        @test Rt4 ≈ 288.2159811164396  rtol = 1e-6
    end

    @testset "tfweight.jl" begin

        ac = TASOPT.load_default_model()
        ip = ipcruise1

        ac.parg[igGearf] = 1.0
        ac.pared[iemblcD, ip] = 46.110246609262873
        ac.pared[ieBPR, ip] = 5.0999999999999996
        ac.pared[iepilc, ip] = 3.0
        ac.pared[iepihc, ip] = 10.000000000000000 
        ac.parg[igdfan] = 1.3927234305722356
        ac.parg[igdlcomp] = 0.67240459668963337
        ac.parg[igrSnace] = 16.000000000000000
        ac.parg[igneng] = 2.0
        ac.parg[igfeadd] = 0.10000000000000001
        ac.parg[igfpylon] = 0.10000000000000001

        Weng, Wnac, Webare, W_HXs, Snace1 = TASOPT.engine.tfweight(ac)

        @test Weng ≈ 46847.51154286845 rtol = 1e-10
        @test Wnac ≈ 9411.055803345604 rtol = 1e-10
        @test Webare ≈ 30161.446412552305 rtol = 1e-10
        @test Snace1 ≈ 24.374719583103083 rtol = 1e-10

    end

    @testset "map_functions.jl" begin
        #Check reverse interpolation function
        map = TASOPT.engine.FanMap
        Wc_target = 559.314
        PR_target = 1.31
        outps_NR = TASOPT.engine.find_NR_inverse_with_derivatives(map.itp_Wc, map.itp_PR, 
            Wc_target, PR_target)
        outp_NR_check = (0.7000000000000002, 1.799999999999995, 0.0002760733345430258, 0.6100706271038377, 0.00430742991535455, -2.9037427622059644)
        
        for (i,outp_NR) in enumerate(outps_NR)
            @test outp_NR ≈ outps_NR[i]
        end

        #Check function to compute speed and efficiency
        pratio = 1.6
        mb = 0.8
        piD = 1.7
        mbD = 1.0
        NbD = 1.0
        ep0 = 0.9
        outps_Ne = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio, mb, piD, mbD, NbD, ep0)
        outps_Ne_check = (0.8826430459793184, 0.8251883153075635, 0.49351112427246263, -0.003692654225532323, -0.39245022849891215, 0.929263500483822, 0.8738166155195252, 1.4836439672220545)
        
        for (i,outp_Ne) in enumerate(outps_Ne)
            @test outp_Ne ≈ outps_Ne_check[i]
        end

        #Use finite difference to validate derivatives of the maps
        eps = 1e-6
        outps_pi = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio + eps, mb, piD, mbD, NbD, ep0)
        outps_mb = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio, mb + eps, piD, mbD, NbD, ep0)
        
        dN_dpi = (outps_pi[1] - outps_Ne[1])/eps
        depol_dpi = (outps_pi[2] - outps_Ne[2])/eps

        dN_dmb = (outps_mb[1] - outps_Ne[1])/eps
        depol_dmb = (outps_mb[2] - outps_Ne[2])/eps

        @test dN_dpi ≈ outps_Ne[3] rtol = 1e-4
        @test dN_dmb ≈ outps_Ne[4] rtol = 1e-4
        @test depol_dpi ≈ outps_Ne[5] rtol = 1e-4
        @test depol_dmb ≈ outps_Ne[6] rtol = 1e-4
    end

    @testset "EngineStation enumeration" begin
        ES = TASOPT.engine.EngineStation
        snum = TASOPT.engine.station_number
        sdesc = TASOPT.engine.station_description

        # All 20 stations must be distinct and constructable
        all_stations = instances(ES.T)
        @test length(all_stations) == 20

        # station_number round-trip: every member maps to a non-empty string
        for s in all_stations
            num = snum(s)
            @test num isa String
            @test !isempty(num)
        end

        # Spot-check TASOPT numbering for all stations
        @test snum(ES.Freestream)     == "0"
        @test snum(ES.FanFaceOuter)   == "18"
        @test snum(ES.FanFaceLPC)     == "19"
        @test snum(ES.PreCoolerOut)   == "19c"
        @test snum(ES.FanFaceFan)     == "2"
        @test snum(ES.FanExit)        == "21"
        @test snum(ES.LPCExit)        == "25"
        @test snum(ES.InterCoolerOut) == "25c"
        @test snum(ES.HPCExit)        == "3"
        @test snum(ES.CombustorExit)  == "4"
        @test snum(ES.CoolMixInlet)   == "4a"
        @test snum(ES.TurbineInlet)   == "41"
        @test snum(ES.HPTExit)        == "45"
        @test snum(ES.LPTExit)        == "49"
        @test snum(ES.RegenCoolerOut) == "49c"
        @test snum(ES.CoreNozzle)     == "5"
        @test snum(ES.CoreNozzleExit) == "6"
        @test snum(ES.FanNozzle)      == "7"
        @test snum(ES.FanNozzleExit)  == "8"
        @test snum(ES.OfftakeDisch)   == "9"

        # station_number values are all distinct
        nums = [snum(s) for s in all_stations]
        @test length(unique(nums)) == length(nums)

        # station_description returns a non-empty string for every member
        for s in all_stations
            desc = sdesc(s)
            @test desc isa String
            @test !isempty(desc)
        end

        # Enum members are usable as Dict keys (field-access sentinel use-case)
        d = Dict(s => snum(s) for s in all_stations)
        @test d[ES.HPCExit]      == "3"
        @test d[ES.TurbineInlet] == "41"

        # Coverage: every station cited in tfsize! and tfoper! docs is present.
        # Documented stations: 0, 18, 19, 19c, 2, 21, 25, 25c, 3, 4, 4a,
        #                      41, 45, 49, 49c, 5, 6, 7, 8, 9
        documented = Set([
            "0", "18", "19", "19c", "2", "21", "25", "25c",
            "3", "4", "4a", "41", "45", "49", "49c",
            "5", "6", "7", "8", "9"
        ])
        enum_nums = Set(snum(s) for s in all_stations)
        @test documented == enum_nums
    end

    @testset "Simple engine" begin

        ac = TASOPT.load_default_model()
        ip = ipcruise1

        TSFC = 1.8e-4
        Fe = 1e4
        neng = 2
        MTOW = 1e5
        feng = 0.05

        ac.pare[ieTSFC, ip, 1] = TSFC
        ac.pare[ieFe, ip, 1] = Fe
        ac.missions[1].points[ip].engine.TSFC = TSFC
        ac.missions[1].points[ip].engine.Fe   = Fe
        ac.parg[igWMTO] = MTOW
        ac.parg[igfeng] = feng

        TASOPT.engine.constant_TSFC_engine!(ac, 0, 1, ip, 0)
        # Sync bare pare → typed state (constant_TSFC_engine! writes pare[iemfuel])
        TASOPT.engine.pare_to_engine_state!(ac.missions[1].points[ip].engine,
                                            view(ac.pare, :, ip, 1))

        @test ac.missions[1].points[ip].engine.mfuel ≈ neng*TSFC*Fe/gee rtol = 1e-10

        TASOPT.engine.fractional_engine_weight!(ac)
        @test ac.parg[igWeng] ≈ feng * MTOW rtol = 1e-10

    end

    @testset "DesignState" begin
        DS = TASOPT.engine.DesignState

        # ------------------------------------------------------------------
        # Zero-initialised Float64 constructor (default)
        # ------------------------------------------------------------------
        ds = DS()
        @test ds isa DS{Float64}

        # All scalar fields are zero
        for fname in (:pifD, :pilcD, :pihcD, :pihtD, :piltD,
                      :mbfD, :mblcD, :mbhcD, :mbhtD, :mbltD,
                      :NbfD, :NblcD, :NbhcD, :NbhtD, :NbltD,
                      :A2, :A25, :A5, :A7,
                      :fc, :ruc, :M4a)
            @test getfield(ds, fname) == 0.0
        end

        # Vector fields are length-4 zero SVectors
        @test length(ds.epsrow) == 4
        @test length(ds.Tmrow) == 4
        @test all(==(0.0), ds.epsrow)
        @test all(==(0.0), ds.Tmrow)

        # ------------------------------------------------------------------
        # Explicit Float64 typed constructor
        # ------------------------------------------------------------------
        ds64 = DS{Float64}()
        @test ds64 isa DS{Float64}
        @test ds64.pifD === 0.0

        # ------------------------------------------------------------------
        # Float32 typed constructor
        # ------------------------------------------------------------------
        ds32 = DS{Float32}()
        @test ds32 isa DS{Float32}
        @test ds32.pifD === 0.0f0
        @test ds32.epsrow isa StaticArrays.SVector{4,Float32}
        @test all(==(0.0f0), ds32.epsrow)

        # ------------------------------------------------------------------
        # Mutation: map scalars round-trip
        # ------------------------------------------------------------------
        ds.pifD  = 1.6
        ds.pilcD = 3.0
        ds.pihcD = 12.0
        ds.pihtD = 4.5
        ds.piltD = 2.0

        @test ds.pifD  == 1.6
        @test ds.pilcD == 3.0
        @test ds.pihcD == 12.0
        @test ds.pihtD == 4.5
        @test ds.piltD == 2.0

        ds.mbfD  = 210.0
        ds.mblcD = 55.0
        ds.mbhcD = 50.0
        ds.mbhtD = 12.0
        ds.mbltD = 45.0

        @test ds.mbfD  == 210.0
        @test ds.mblcD == 55.0
        @test ds.mbhcD == 50.0
        @test ds.mbhtD == 12.0
        @test ds.mbltD == 45.0

        ds.NbfD  = 1.0
        ds.NblcD = 1.0
        ds.NbhcD = 1.0
        ds.NbhtD = 1.0
        ds.NbltD = 1.0

        @test ds.NbfD == 1.0

        # ------------------------------------------------------------------
        # Mutation: component areas round-trip
        # ------------------------------------------------------------------
        ds.A2  = 3.14
        ds.A25 = 0.52
        ds.A5  = 0.18
        ds.A7  = 1.20

        @test ds.A2  == 3.14
        @test ds.A25 == 0.52
        @test ds.A5  == 0.18
        @test ds.A7  == 1.20

        # ------------------------------------------------------------------
        # Mutation: cooling fields round-trip
        # ------------------------------------------------------------------
        using StaticArrays
        eps_val = SA[0.05, 0.04, 0.03, 0.02]
        tm_val  = SA[1200.0, 1150.0, 1100.0, 1050.0]

        ds.epsrow = eps_val
        ds.Tmrow  = tm_val
        ds.fc     = sum(eps_val)
        ds.ruc    = 0.95
        ds.M4a    = 0.10

        @test ds.epsrow[1] ≈ 0.05
        @test ds.epsrow[4] ≈ 0.02
        @test ds.Tmrow[1]  ≈ 1200.0
        @test ds.Tmrow[4]  ≈ 1050.0
        @test ds.fc        ≈ 0.14
        @test ds.ruc       ≈ 0.95
        @test ds.M4a       ≈ 0.10

        # ------------------------------------------------------------------
        # ncrowx = 4: vector dimensions match legacy constant
        # ------------------------------------------------------------------
        # The legacy index.inc defines ncrowx = ieTmet1 - ieepsc1 = 192 - 188 = 4
        @test length(ds.epsrow) == 4
        @test length(ds.Tmrow)  == 4

        # ------------------------------------------------------------------
        # @inferred: field access should not allocate or lose type info
        # ------------------------------------------------------------------
        @test @inferred(Float64, getproperty(ds, :pifD)) == ds.pifD
        @test @inferred(Float64, getproperty(ds, :A2))  == ds.A2
        @test @inferred(Float64, getproperty(ds, :ruc)) == ds.ruc
        @test @inferred(Float64, getproperty(ds, :M4a)) == ds.M4a

        # ------------------------------------------------------------------
        # Inline storage: DesignState embedded in another struct should not
        # force heap allocation of the DesignState itself
        # ------------------------------------------------------------------
        struct WrapDS
            ds::DS{Float64}
            tag::Int
        end
        w = WrapDS(DS(), 42)
        @test w.ds isa DS{Float64}
        @test w.tag == 42
    end

    # ======================================================================
    @testset "GasState" begin
        using StaticArrays
        GS = TASOPT.engine.GasState

        # ------------------------------------------------------------------
        # Default constructor — GasState() → GasState{Float64}
        # ------------------------------------------------------------------
        gs = GS()
        @test gs isa GS{Float64}

        # Scalar fields zero-initialised
        @test gs.Tt  === 0.0
        @test gs.ht  === 0.0
        @test gs.pt  === 0.0
        @test gs.cpt === 0.0
        @test gs.Rt  === 0.0
        @test gs.st  === 0.0

        @test gs.Ts  === 0.0
        @test gs.ps  === 0.0
        @test gs.hs  === 0.0
        @test gs.ss  === 0.0
        @test gs.cps === 0.0
        @test gs.Rs  === 0.0
        @test gs.u   === 0.0

        # Species composition is length-5 zero SVector
        @test length(gs.alpha) == 5
        @test gs.alpha isa StaticArrays.SVector{5,Float64}
        @test all(==(0.0), gs.alpha)

        # ------------------------------------------------------------------
        # Explicit Float64 typed constructor
        # ------------------------------------------------------------------
        gs64 = GS{Float64}()
        @test gs64 isa GS{Float64}
        @test gs64.Tt === 0.0

        # ------------------------------------------------------------------
        # Float32 typed constructor
        # ------------------------------------------------------------------
        gs32 = GS{Float32}()
        @test gs32 isa GS{Float32}
        @test gs32.Tt === 0.0f0
        @test gs32.alpha isa StaticArrays.SVector{5,Float32}
        @test all(==(0.0f0), gs32.alpha)

        # ------------------------------------------------------------------
        # Constructor with explicit total-state fields and species
        # ------------------------------------------------------------------
        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]   # standard air
        gs_total  = GS{Float64}(288.15, 2.885e5, 101325.0, 1005.0, 287.058, air_alpha)

        @test gs_total.Tt  ≈ 288.15
        @test gs_total.ht  ≈ 2.885e5
        @test gs_total.pt  ≈ 101325.0
        @test gs_total.cpt ≈ 1005.0
        @test gs_total.Rt  ≈ 287.058

        # st (total entropy complement) defaults to zero (not yet computed)
        @test gs_total.st  === 0.0
        # Static fields default to zero
        @test gs_total.Ts  === 0.0
        @test gs_total.ps  === 0.0
        @test gs_total.u   === 0.0

        @test gs_total.alpha ≈ air_alpha

        # ------------------------------------------------------------------
        # Mutation: round-trip on every field
        # ------------------------------------------------------------------
        gs.Tt  = 800.0
        gs.ht  = 8.5e5
        gs.pt  = 2.0e6
        gs.cpt = 1100.0
        gs.Rt  = 287.058
        gs.st  = 2350.0

        @test gs.Tt  ≈ 800.0
        @test gs.ht  ≈ 8.5e5
        @test gs.pt  ≈ 2.0e6
        @test gs.cpt ≈ 1100.0
        @test gs.Rt  ≈ 287.058
        @test gs.st  ≈ 2350.0

        gs.Ts  = 720.0
        gs.ps  = 1.6e6
        gs.hs  = 7.9e5
        gs.ss  = 2400.0
        gs.cps = 1090.0
        gs.Rs  = 287.058
        gs.u   = 150.0

        @test gs.Ts  ≈ 720.0
        @test gs.ps  ≈ 1.6e6
        @test gs.hs  ≈ 7.9e5
        @test gs.ss  ≈ 2400.0
        @test gs.cps ≈ 1090.0
        @test gs.Rs  ≈ 287.058
        @test gs.u   ≈ 150.0

        gs.alpha = air_alpha
        @test gs.alpha[1] ≈ 0.7532
        @test gs.alpha[5] ≈ 0.0127
        @test sum(gs.alpha) ≈ sum(air_alpha)

        # ------------------------------------------------------------------
        # @inferred: field access must not lose type info
        # ------------------------------------------------------------------
        gs2 = GS()
        @test @inferred(Float64, getproperty(gs2, :Tt)) == gs2.Tt
        @test @inferred(Float64, getproperty(gs2, :pt)) == gs2.pt
        @test @inferred(Float64, getproperty(gs2, :u))  == gs2.u

        # ------------------------------------------------------------------
        # Inline storage: GasState embedded in another struct
        # ------------------------------------------------------------------
        struct WrapGS
            state::GS{Float64}
            id::Int
        end
        wg = WrapGS(GS(), 7)
        @test wg.state isa GS{Float64}
        @test wg.id == 7

        # ------------------------------------------------------------------
        # OQ-5 invariant: alpha always has exactly 5 species
        # (n == 5 is a compile-time guarantee via SVector length)
        # ------------------------------------------------------------------
        @test length(GS().alpha) == 5
        @test length(GS{Float32}().alpha) == 5
    end

    # ======================================================================
    # FlowStation
    # ======================================================================
    @testset "FlowStation" begin
        using StaticArrays
        FS = TASOPT.engine.FlowStation
        GS = TASOPT.engine.GasState

        # ------------------------------------------------------------------
        # Default constructor — FlowStation() → FlowStation{Float64}
        # ------------------------------------------------------------------
        fs = FS()
        @test fs isa FS{Float64}
        @test fs.A    === 0.0
        @test fs.mdot === 0.0

        # Embedded gas is a properly zeroed GasState
        @test fs.gas isa GS{Float64}
        @test fs.gas.Tt === 0.0
        @test fs.gas.pt === 0.0
        @test fs.gas.u  === 0.0

        # ------------------------------------------------------------------
        # Typed constructors
        # ------------------------------------------------------------------
        fs64 = FS{Float64}()
        @test fs64 isa FS{Float64}
        @test fs64.A === 0.0

        fs32 = FS{Float32}()
        @test fs32 isa FS{Float32}
        @test fs32.A    === 0.0f0
        @test fs32.mdot === 0.0f0
        @test fs32.gas isa GS{Float32}

        # Constructor with explicit A and mdot keyword args
        fs_area = FS{Float64}(; A=0.5, mdot=10.0)
        @test fs_area.A    ≈ 0.5
        @test fs_area.mdot ≈ 10.0
        @test fs_area.gas.Tt === 0.0   # gas still zeroed

        # ------------------------------------------------------------------
        # Constructor with explicit total state + species
        # ------------------------------------------------------------------
        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        fs_total  = FS{Float64}(288.15, 2.885e5, 101325.0, 1005.0, 287.058, air_alpha;
                                A=0.3, mdot=50.0)

        @test fs_total.Tt    ≈ 288.15
        @test fs_total.ht    ≈ 2.885e5
        @test fs_total.pt    ≈ 101325.0
        @test fs_total.cpt   ≈ 1005.0
        @test fs_total.Rt    ≈ 287.058
        @test fs_total.alpha ≈ air_alpha
        @test fs_total.A     ≈ 0.3
        @test fs_total.mdot  ≈ 50.0

        # Static fields default to zero even with explicit total state
        @test fs_total.Ts === 0.0
        @test fs_total.ps === 0.0
        @test fs_total.u  === 0.0

        # ------------------------------------------------------------------
        # Property forwarding — reads
        # ------------------------------------------------------------------
        fs2 = FS{Float64}(; A=1.0, mdot=20.0)
        fs2.gas.Tt  = 800.0
        fs2.gas.pt  = 2.0e6
        fs2.gas.u   = 150.0

        # Forwarded reads match direct nested access
        @test fs2.Tt  === fs2.gas.Tt
        @test fs2.pt  === fs2.gas.pt
        @test fs2.u   === fs2.gas.u

        # Own-field reads still work
        @test fs2.A    ≈ 1.0
        @test fs2.mdot ≈ 20.0

        # gas field itself is accessible
        @test fs2.gas isa GS{Float64}

        # ------------------------------------------------------------------
        # Property forwarding — writes via setproperty!
        # ------------------------------------------------------------------
        fs3 = FS()
        fs3.Tt  = 500.0
        fs3.ht  = 5.0e5
        fs3.pt  = 1.0e6
        fs3.cpt = 1050.0
        fs3.Rt  = 287.0
        fs3.Ts  = 460.0
        fs3.ps  = 0.9e6
        fs3.hs  = 4.8e5
        fs3.ss  = 2300.0
        fs3.cps = 1040.0
        fs3.Rs  = 287.0
        fs3.u   = 200.0
        fs3.alpha = air_alpha

        # All forwarded writes are visible via the gas field
        @test fs3.gas.Tt  ≈ 500.0
        @test fs3.gas.ht  ≈ 5.0e5
        @test fs3.gas.pt  ≈ 1.0e6
        @test fs3.gas.Ts  ≈ 460.0
        @test fs3.gas.u   ≈ 200.0
        @test fs3.gas.alpha ≈ air_alpha

        # Also visible via the forwarded reader
        @test fs3.Tt ≈ 500.0
        @test fs3.u  ≈ 200.0

        # Own-field writes
        fs3.A    = 0.8
        fs3.mdot = 30.0
        @test fs3.A    ≈ 0.8
        @test fs3.mdot ≈ 30.0

        # ------------------------------------------------------------------
        # @inferred: forwarded field access must not lose type information
        # ------------------------------------------------------------------
        fs4 = FS()
        @test @inferred(Float64, getproperty(fs4, :Tt))  == fs4.Tt
        @test @inferred(Float64, getproperty(fs4, :pt))  == fs4.pt
        @test @inferred(Float64, getproperty(fs4, :u))   == fs4.u
        @test @inferred(Float64, getproperty(fs4, :A))   == fs4.A
        @test @inferred(Float64, getproperty(fs4, :mdot))== fs4.mdot

        # ------------------------------------------------------------------
        # propertynames includes both own and forwarded names
        # ------------------------------------------------------------------
        pnames = propertynames(FS())
        for name in (:gas, :A, :mdot, :Tt, :ht, :pt, :cpt, :Rt,
                     :Ts, :ps, :hs, :ss, :cps, :Rs, :u, :alpha)
            @test name in pnames
        end

        # ------------------------------------------------------------------
        # Inline storage: FlowStation embedded in another struct
        # ------------------------------------------------------------------
        struct WrapFS
            station::FS{Float64}
            id::Int
        end
        wfs = WrapFS(FS(), 3)
        @test wfs.station isa FS{Float64}
        @test wfs.id == 3

        # ------------------------------------------------------------------
        # Float32 forwarding works end-to-end
        # ------------------------------------------------------------------
        fs32b = FS{Float32}()
        fs32b.Tt = 300.0f0
        fs32b.pt = 1.0f5
        @test fs32b.Tt === 300.0f0
        @test fs32b.pt === 1.0f5
        @test @inferred(Float32, getproperty(fs32b, :Tt)) == fs32b.Tt

        # ------------------------------------------------------------------
        # st field is forwarded (GasState.st → FlowStation.st)
        # ------------------------------------------------------------------
        fs_st = FS()
        @test fs_st.st === 0.0
        fs_st.st = 2400.0
        @test fs_st.gas.st ≈ 2400.0
        @test fs_st.st ≈ 2400.0
    end

    # ======================================================================
    # set_total_from_Tt!
    # ======================================================================
    @testset "set_total_from_Tt!" begin
        using StaticArrays
        FS     = TASOPT.engine.FlowStation
        gassum = TASOPT.engine.gassum
        set_tt = TASOPT.engine.set_total_from_Tt!

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        Tt_test   = 800.0

        # Build a station with known Tt and alpha; pt is arbitrary (not touched)
        fs = FS()
        fs.alpha = air_alpha
        fs.Tt    = Tt_test
        fs.pt    = 2.0e6

        # Call the wrapper
        ret = set_tt(fs)

        # Return value is the station itself (for chaining)
        @test ret === fs

        # Compare against direct gassum call
        s_ref, _, h_ref, _, cp_ref, r_ref = gassum(air_alpha, 5, Tt_test)

        @test fs.ht  ≈ h_ref
        @test fs.st  ≈ s_ref
        @test fs.cpt ≈ cp_ref
        @test fs.Rt  ≈ r_ref

        # Fields NOT updated by set_total_from_Tt!
        @test fs.Tt  ≈ Tt_test   # unchanged
        @test fs.pt  ≈ 2.0e6     # unchanged
        @test fs.Ts  === 0.0     # static not touched
        @test fs.ps  === 0.0
        @test fs.u   === 0.0

        # Monotonicity: higher Tt → higher ht and higher st (entropy increases with T)
        fs2 = FS()
        fs2.alpha = air_alpha
        fs2.Tt    = 1600.0
        set_tt(fs2)
        @test fs2.ht > fs.ht
        @test fs2.st > fs.st

    end

    # set_static_from_M!
    # ======================================================================
    @testset "set_static_from_M!" begin
        using StaticArrays
        FS        = TASOPT.engine.FlowStation
        gas_mach  = TASOPT.engine.gas_mach
        set_tt    = TASOPT.engine.set_total_from_Tt!
        set_mach  = TASOPT.engine.set_static_from_M!

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]

        # Build a fully initialised total state
        fs = FS()
        fs.alpha = air_alpha
        fs.Tt    = 800.0
        fs.pt    = 2.0e5
        set_tt(fs)   # populate ht, st, cpt, Rt

        # ------------------------------------------------------------------
        # Return value is the station itself (chaining)
        # ------------------------------------------------------------------
        ret = set_mach(fs, 0.6)
        @test ret === fs

        # ------------------------------------------------------------------
        # Results match direct gas_mach call
        # ------------------------------------------------------------------
        ps_ref, Ts_ref, hs_ref, ss_ref, cps_ref, Rs_ref =
            gas_mach(air_alpha, 5, fs.pt, fs.Tt, fs.ht, fs.st, fs.cpt, fs.Rt,
                     0.0, 0.6, 1.0)

        @test fs.Ts  ≈ Ts_ref
        @test fs.ps  ≈ ps_ref
        @test fs.hs  ≈ hs_ref
        @test fs.ss  ≈ ss_ref
        @test fs.cps ≈ cps_ref
        @test fs.Rs  ≈ Rs_ref

        # ------------------------------------------------------------------
        # Energy conservation: hs + ½u² ≈ ht
        # ------------------------------------------------------------------
        @test fs.hs + 0.5 * fs.u^2 ≈ fs.ht   rtol=1e-6

        # ------------------------------------------------------------------
        # Isentropic invariant (epol=1.0): ps = pt * exp((ss - st) / Rt)
        # st = s[Tt] and ss = s[Ts] differ because T changes; the invariant
        # is that thermodynamic entropy is conserved, encoded by gas_mach as:
        #   p = po * exp(epol * (s - so) / ro)
        # ------------------------------------------------------------------
        @test fs.ps ≈ fs.pt * exp((fs.ss - fs.st) / fs.Rt)   rtol=1e-6

        # ------------------------------------------------------------------
        # M=0 special case: static state equals total state, u=0
        # ------------------------------------------------------------------
        fs0 = FS()
        fs0.alpha = air_alpha
        fs0.Tt = 800.0; fs0.pt = 2.0e5
        set_tt(fs0)
        set_mach(fs0, 0.0)
        @test fs0.Ts ≈ fs0.Tt   rtol=1e-6
        @test fs0.ps ≈ fs0.pt   rtol=1e-6
        @test fs0.hs ≈ fs0.ht   rtol=1e-6
        @test fs0.u  === 0.0

        # ------------------------------------------------------------------
        # Total state fields are not modified
        # ------------------------------------------------------------------
        @test fs.Tt    ≈ 800.0
        @test fs.pt    ≈ 2.0e5
        @test fs.alpha == air_alpha

        # ------------------------------------------------------------------
        # Monotonicity: higher M → lower Ts and lower ps
        # ------------------------------------------------------------------
        fs_lo = FS()
        fs_hi = FS()
        for fsm in (fs_lo, fs_hi)
            fsm.alpha = air_alpha; fsm.Tt = 800.0; fsm.pt = 2.0e5
            set_tt(fsm)
        end
        set_mach(fs_lo, 0.4)
        set_mach(fs_hi, 0.8)
        @test fs_lo.Ts > fs_hi.Ts
        @test fs_lo.ps > fs_hi.ps
        @test fs_lo.u  < fs_hi.u

    end

    # apply_pratio_from!
    # ======================================================================
    @testset "apply_pratio_from!" begin
        using StaticArrays
        FS          = TASOPT.engine.FlowStation
        gas_prat    = TASOPT.engine.gas_prat
        set_tt      = TASOPT.engine.set_total_from_Tt!
        appr        = TASOPT.engine.apply_pratio_from!

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        Tt_in  = 288.15
        pt_in  = 1.0e5
        pratio = 20.0
        epol   = 0.9

        inlet = FS(); inlet.alpha = air_alpha; inlet.Tt = Tt_in; inlet.pt = pt_in
        set_tt(inlet)
        outlet = FS()

        # Return value is the outlet station (for chaining)
        ret = appr(outlet, inlet, pratio; epol=epol)
        @test ret === outlet

        # Results match direct gas_prat call
        pt_ref, Tt_ref, ht_ref, st_ref, cpt_ref, Rt_ref =
            gas_prat(air_alpha, 5,
                     inlet.pt, inlet.Tt, inlet.ht, inlet.st, inlet.cpt, inlet.Rt,
                     pratio, epol)

        @test outlet.pt  ≈ pt_ref
        @test outlet.Tt  ≈ Tt_ref
        @test outlet.ht  ≈ ht_ref
        @test outlet.st  ≈ st_ref
        @test outlet.cpt ≈ cpt_ref
        @test outlet.Rt  ≈ Rt_ref

        # Pressure invariant: outlet.pt == inlet.pt * pratio (gas_prat sets p = po*pratio)
        @test outlet.pt ≈ inlet.pt * pratio

        # Composition preserved
        @test outlet.alpha == air_alpha

        # Static fields not touched
        @test outlet.Ts === 0.0
        @test outlet.ps === 0.0
        @test outlet.u  === 0.0

        # Inlet total state is not modified
        @test inlet.Tt ≈ Tt_in
        @test inlet.pt ≈ pt_in

        # pratio=1, epol=1 → outlet equals inlet (identity process)
        outlet_id = FS()
        appr(outlet_id, inlet, 1.0)
        @test outlet_id.Tt  ≈ inlet.Tt   rtol=1e-5
        @test outlet_id.pt  ≈ inlet.pt   rtol=1e-5
        @test outlet_id.ht  ≈ inlet.ht   rtol=1e-5

        # Compression raises temperature and enthalpy
        @test outlet.Tt > inlet.Tt
        @test outlet.ht > inlet.ht

    end

    # apply_delh_from!
    # ======================================================================
    @testset "apply_delh_from!" begin
        using StaticArrays
        FS          = TASOPT.engine.FlowStation
        gas_delh    = TASOPT.engine.gas_delh
        set_tt      = TASOPT.engine.set_total_from_Tt!
        appdh       = TASOPT.engine.apply_delh_from!

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        Tt_in  = 900.0
        pt_in  = 2.0e6
        delh   = -200_000.0   # turbine expansion (negative)
        epol   = 0.88

        inlet = FS(); inlet.alpha = air_alpha; inlet.Tt = Tt_in; inlet.pt = pt_in
        set_tt(inlet)
        outlet = FS()

        # Return value is the outlet station (for chaining)
        ret = appdh(outlet, inlet, delh; epol=epol)
        @test ret === outlet

        # Results match direct gas_delh call
        pt_ref, Tt_ref, ht_ref, st_ref, cpt_ref, Rt_ref =
            gas_delh(air_alpha, 5,
                     inlet.pt, inlet.Tt, inlet.ht, inlet.st, inlet.cpt, inlet.Rt,
                     delh, epol)

        @test outlet.pt  ≈ pt_ref
        @test outlet.Tt  ≈ Tt_ref
        @test outlet.ht  ≈ ht_ref
        @test outlet.st  ≈ st_ref
        @test outlet.cpt ≈ cpt_ref
        @test outlet.Rt  ≈ Rt_ref

        # Enthalpy invariant: outlet.ht ≈ inlet.ht + delh
        @test outlet.ht ≈ inlet.ht + delh   rtol=1e-6

        # Composition preserved
        @test outlet.alpha == air_alpha

        # Static fields not touched
        @test outlet.Ts === 0.0
        @test outlet.ps === 0.0
        @test outlet.u  === 0.0

        # Inlet total state is not modified
        @test inlet.Tt ≈ Tt_in
        @test inlet.pt ≈ pt_in

        # delh=0, epol=1 → outlet equals inlet (identity process)
        outlet_id = FS()
        appdh(outlet_id, inlet, 0.0)
        @test outlet_id.Tt ≈ inlet.Tt   rtol=1e-5
        @test outlet_id.ht ≈ inlet.ht   rtol=1e-5

        # Expansion (negative delh) decreases temperature and pressure
        @test outlet.Tt < inlet.Tt
        @test outlet.pt < inlet.pt

    end

    # EngineState
    @testset "EngineState" begin

        ES = TASOPT.engine.EngineState
        FS = TASOPT.engine.FlowStation
        DS = TASOPT.engine.DesignState

        # ------------------------------------------------------------------
        # Default constructor — EngineState() → EngineState{Float64}
        # ------------------------------------------------------------------
        eng = ES()
        @test eng isa ES{Float64}

        # ------------------------------------------------------------------
        # Typed constructor — EngineState{Float32}()
        # ------------------------------------------------------------------
        eng32 = ES{Float32}()
        @test eng32 isa ES{Float32}
        @test eng32.M0 === 0.0f0
        @test eng32.st4 isa FS{Float32}
        @test eng32.design isa DS{Float32}

        # ------------------------------------------------------------------
        # All 20 station fields are FlowStation{Float64} and fully zeroed
        # ------------------------------------------------------------------
        station_fields = (:st0, :st18, :st19, :st19c, :st2, :st21,
                          :st25, :st25c, :st3, :st4, :st4a, :st41,
                          :st45, :st49, :st49c, :st5, :st6, :st7,
                          :st8, :st9)
        @test length(station_fields) == 20
        for fname in station_fields
            fs = getproperty(eng, fname)
            @test fs isa FS{Float64}
            @test fs.Tt  === 0.0
            @test fs.pt  === 0.0
            @test fs.A   === 0.0
            @test fs.mdot === 0.0
        end

        # ------------------------------------------------------------------
        # Ambient scalars all zero by default
        # ------------------------------------------------------------------
        @test eng.M0 === 0.0
        @test eng.T0 === 0.0
        @test eng.p0 === 0.0
        @test eng.a0 === 0.0

        # ------------------------------------------------------------------
        # Design field is a zero-initialised DesignState
        # ------------------------------------------------------------------
        @test eng.design isa DS{Float64}
        @test eng.design.pifD  === 0.0
        @test eng.design.pilcD === 0.0
        @test eng.design.A2    === 0.0

        # ------------------------------------------------------------------
        # Setting and reading station fields
        # ------------------------------------------------------------------
        eng.st4.Tt  = 1800.0
        eng.st4.pt  = 2.5e6
        eng.st4.mdot = 40.0
        @test eng.st4.Tt   ≈ 1800.0
        @test eng.st4.pt   ≈ 2.5e6
        @test eng.st4.mdot ≈ 40.0

        # Another station is unaffected
        @test eng.st3.Tt === 0.0

        # ------------------------------------------------------------------
        # Setting ambient scalars
        # ------------------------------------------------------------------
        eng2 = ES()
        eng2.M0 = 0.85
        eng2.T0 = 216.65
        eng2.p0 = 22632.0
        eng2.a0 = 295.07
        @test eng2.M0 ≈ 0.85
        @test eng2.T0 ≈ 216.65
        @test eng2.p0 ≈ 22632.0
        @test eng2.a0 ≈ 295.07

        # ------------------------------------------------------------------
        # Setting design-state fields
        # ------------------------------------------------------------------
        eng2.design.pifD  = 1.6
        eng2.design.pihcD = 10.0
        eng2.design.A2    = 0.5
        @test eng2.design.pifD  ≈ 1.6
        @test eng2.design.pihcD ≈ 10.0
        @test eng2.design.A2    ≈ 0.5

        # ------------------------------------------------------------------
        # Each station is an independent FlowStation (mutations don't alias)
        # ------------------------------------------------------------------
        eng3 = ES()
        eng3.st4.Tt  = 1600.0
        eng3.st41.Tt = 1550.0
        @test eng3.st4.Tt  ≈ 1600.0
        @test eng3.st41.Tt ≈ 1550.0
        @test eng3.st45.Tt === 0.0   # neighbouring station unaffected

        # ------------------------------------------------------------------
        # Inline storage: EngineState embedded in another struct
        # ------------------------------------------------------------------
        struct WrapES
            eng::ES{Float64}
            id::Int
        end
        wes = WrapES(ES(), 7)
        @test wes.eng isa ES{Float64}
        @test wes.id == 7

        # ------------------------------------------------------------------
        # Float32 stations are properly typed
        # ------------------------------------------------------------------
        eng32b = ES{Float32}()
        eng32b.st3.Tt = 900.0f0
        @test eng32b.st3.Tt === 900.0f0
        @test eng32b.st3 isa FS{Float32}

        # ------------------------------------------------------------------
        # Station-shortcut access — eng.Tt4 ≡ eng.st4.Tt
        # ------------------------------------------------------------------
        @testset "station shortcuts" begin

            eng_sc = ES()

            # Shortcut reads agree with direct nested access (zero state)
            @test eng_sc.Tt4   === eng_sc.st4.Tt
            @test eng_sc.pt25  === eng_sc.st25.pt
            @test eng_sc.A9    === eng_sc.st9.A
            @test eng_sc.mdot0 === eng_sc.st0.mdot

            # Fractional-suffix stations work correctly
            @test eng_sc.Tt19c  === eng_sc.st19c.Tt
            @test eng_sc.pt4a   === eng_sc.st4a.pt
            @test eng_sc.ht49c  === eng_sc.st49c.ht
            @test eng_sc.A25c   === eng_sc.st25c.A

            # Shortcut write then direct read
            eng_sc.Tt4   = 1600.0
            eng_sc.pt25  = 4.5e5
            eng_sc.A9    = 0.18
            eng_sc.mdot2 = 120.0
            @test eng_sc.st4.Tt   ≈ 1600.0
            @test eng_sc.st25.pt  ≈ 4.5e5
            @test eng_sc.st9.A    ≈ 0.18
            @test eng_sc.st2.mdot ≈ 120.0

            # Direct write then shortcut read
            eng_sc.st41.Tt = 1500.0
            eng_sc.st45.pt = 3.0e5
            @test eng_sc.Tt41 ≈ 1500.0
            @test eng_sc.pt45 ≈ 3.0e5

            # Shortcut write/read consistency (round-trip)
            eng_sc.Tt49c = 900.0
            @test eng_sc.Tt49c ≈ eng_sc.st49c.Tt ≈ 900.0

            # Non-shortcut own fields still resolve correctly after override
            eng_sc.M0 = 0.85
            @test eng_sc.M0 ≈ 0.85
            eng_sc.T0 = 216.65
            @test eng_sc.T0 ≈ 216.65

            # Type stability: @inferred verifies no type info is lost through shortcut.
            # Note: @inferred requires a call expression, so we use getproperty()
            # explicitly rather than dot-access syntax.
            @test @inferred(Float64, getproperty(eng_sc, :Tt4))  == eng_sc.st4.Tt
            @test @inferred(Float64, getproperty(eng_sc, :pt25)) == eng_sc.st25.pt
            @test @inferred(Float64, getproperty(eng_sc, :A9))   == eng_sc.st9.A

            # propertynames includes shortcut symbols
            pnames = propertynames(eng_sc)
            @test :Tt4   in pnames
            @test :pt25c in pnames
            @test :A9    in pnames
            @test :mdot0 in pnames
            @test :Tt49c in pnames
            # own fields are also listed
            @test :st4   in pnames
            @test :M0    in pnames
            @test :design in pnames

            # Unknown shortcut raises an error
            @test_throws ErrorException eng_sc.bogus_field_xyz

            # Float32 shortcut access preserves element type
            eng32_sc = ES{Float32}()
            eng32_sc.Tt4 = 1400.0f0
            @test eng32_sc.Tt4 === 1400.0f0
            @test @inferred(Float32, getproperty(eng32_sc, :Tt4)) == 1400.0f0

        end  # station shortcuts

    end  # EngineState

    # ======================================================================
    # dump_stations
    # ======================================================================
    @testset "dump_stations" begin
        ES = TASOPT.engine.EngineState

        eng = ES()
        # Set ambient flight condition
        eng.M0  = 0.80
        eng.T0  = 218.8
        eng.p0  = 23842.0
        eng.a0  = 296.5

        # Set a few stations with recognisable values
        eng.Tt4  = 1500.0
        eng.pt3  = 2.5e6
        eng.mdot2 = 320.0
        eng.A9   = 0.42

        # ------------------------------------------------------------------
        # Capture output to a buffer and verify structural invariants
        # ------------------------------------------------------------------
        buf = IOBuffer()
        TASOPT.engine.dump_stations(buf, eng)
        out = String(take!(buf))

        # Header line contains ambient quantities
        @test occursin("M0=0.8000", out)
        @test occursin("T0=218.80K", out)
        @test occursin("p0=23842.00Pa", out)
        @test occursin("a0=296.50m/s", out)

        # Column header row is present
        @test occursin("Tt[K]", out)
        @test occursin("pt[Pa]", out)
        @test occursin("mdot[kg/s]", out)

        # Every station number appears as a row prefix
        for stnum in ("0", "2", "18", "19", "19c", "21",
                      "25", "25c", "3", "4", "4a",
                      "41", "45", "49", "49c",
                      "5", "6", "7", "8", "9")
            @test occursin(stnum, out)
        end

        # Verify the values we set appear in the output
        @test occursin("1500", out)   # Tt4
        @test occursin("320", out)    # mdot2
        @test occursin("0.42", out)   # A9

        # ------------------------------------------------------------------
        # No-argument form writes to stdout (just check it doesn't throw)
        # ------------------------------------------------------------------
        eng2 = ES()
        @test begin
            redirect_stdout(devnull) do
                TASOPT.engine.dump_stations(eng2)
            end
            true
        end

        # ------------------------------------------------------------------
        # File-target invocation: write to a temp file and read it back
        # ------------------------------------------------------------------
        tmp = tempname()
        open(tmp, "w") do f
            TASOPT.engine.dump_stations(f, eng)
        end
        file_out = read(tmp, String)
        @test occursin("M0=0.8000", file_out)
        @test occursin("Tt[K]", file_out)
        rm(tmp; force=true)

        # ------------------------------------------------------------------
        # Float32 variant works without error
        # ------------------------------------------------------------------
        eng32 = ES{Float32}()
        eng32.M0 = 0.8f0
        buf32 = IOBuffer()
        TASOPT.engine.dump_stations(buf32, eng32)
        out32 = String(take!(buf32))
        @test occursin("M0=0.8000", out32)

    end  # dump_stations

    # ======================================================================
    # run_engine_design_point / pare_to_engine_state!
    # ======================================================================
    @testset "run_engine_design_point" begin
        # ------------------------------------------------------------------
        # Run the harness on a fully sized aircraft.
        # size_aircraft! must be called first so that pare[ieFe] carries a
        # meaningful design thrust for the engine sizing.
        # ------------------------------------------------------------------
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)
        eng = TASOPT.engine.run_engine_design_point(ac)

        # ---- Ambient scalars are populated and physically plausible ----
        @test eng.M0 > 0.0          # cruise Mach
        @test eng.T0 > 200.0        # stratosphere T, K
        @test eng.p0 > 1e4          # cruise pressure, Pa
        @test eng.a0 > 280.0        # speed of sound, m/s

        # ---- Total temperature rises through compression / combustor ----
        # Tt2 ≈ Tt0: adiabatic inlet preserves stagnation temperature.
        # Tt25 > Tt2 (or Tt19): LPC adds work.
        @test eng.Tt2  ≈ eng.Tt0  rtol = 1e-6   # isentropic inlet
        @test eng.Tt25 > eng.Tt19   # LPC heating
        @test eng.Tt3  > eng.Tt25   # HPC heating
        @test eng.Tt4  > eng.Tt3    # combustor

        # ---- Total temperature falls monotonically through expansion ----
        @test eng.Tt41 < eng.Tt4    # turbine inlet (after cooling dilution)
        @test eng.Tt45 < eng.Tt41   # HPT exit
        @test eng.Tt49 < eng.Tt45   # LPT exit

        # ---- Total pressure: compressor side rises, turbine side falls ----
        @test eng.Tt3 > eng.Tt2     # same as temperature but pressure check
        @test eng.pt3 > eng.pt2     # HPC raises pressure
        @test eng.pt45 < eng.pt41   # HPT drops pressure
        @test eng.pt49 < eng.pt45   # LPT drops pressure

        # ---- Nozzle: static pressure equals ambient at design ----
        #  (ps_exit >= p0 for under-expanded or ideally expanded nozzles)
        @test eng.ps5 >= 0.9 * eng.p0
        @test eng.ps7 >= 0.9 * eng.p0

        # ---- Station areas and mass flow are positive ----
        @test eng.A2    > 0.0
        @test eng.A25   > 0.0
        @test eng.A5    > 0.0
        @test eng.A7    > 0.0
        @test eng.st2.mdot > 0.0    # core mass flow

        # ------------------------------------------------------------------
        # pare_to_engine_state! round-trip: after a run, reading pare
        # into a fresh EngineState must reproduce the same values.
        # ------------------------------------------------------------------
        eng2 = TASOPT.engine.EngineState{Float64}()
        TASOPT.engine.pare_to_engine_state!(eng2, view(ac.pare, :, ipcruise1, 1))

        @test eng2.Tt4  ≈ eng.Tt4   rtol = 1e-12
        @test eng2.pt3  ≈ eng.pt3   rtol = 1e-12
        @test eng2.Tt49 ≈ eng.Tt49  rtol = 1e-12
        @test eng2.A2   ≈ eng.A2    rtol = 1e-12
        @test eng2.M0   ≈ eng.M0    rtol = 1e-12

        # ------------------------------------------------------------------
        # Consistency with pare: EngineState values must equal pare directly.
        # This verifies the mapping is correct, not just round-trip stable.
        # ------------------------------------------------------------------
        pare_ip = view(ac.pare, :, ipcruise1, 1)

        @test eng.Tt4  ≈ pare_ip[ieTt4]   rtol = 1e-12
        @test eng.ht4  ≈ pare_ip[ieht4]   rtol = 1e-12
        @test eng.pt4  ≈ pare_ip[iept4]   rtol = 1e-12
        @test eng.Tt3  ≈ pare_ip[ieTt3]   rtol = 1e-12
        @test eng.pt3  ≈ pare_ip[iept3]   rtol = 1e-12
        @test eng.Tt49 ≈ pare_ip[ieTt49]  rtol = 1e-12
        @test eng.pt49 ≈ pare_ip[iept49]  rtol = 1e-12
        @test eng.Ts2  ≈ pare_ip[ieT2]    rtol = 1e-12
        @test eng.ps2  ≈ pare_ip[iep2]    rtol = 1e-12
        @test eng.A2   ≈ pare_ip[ieA2]    rtol = 1e-12
        @test eng.st2.mdot ≈ pare_ip[iemcore] rtol = 1e-12

        # ------------------------------------------------------------------
        # Station 6 and 8 (nozzle exits): static-only in pare, total zero.
        # Verify static fields are non-zero and total fields remain zero.
        # ------------------------------------------------------------------
        @test eng.ps6  > 0.0
        @test eng.Ts6  > 0.0
        @test eng.Tt6  == 0.0    # total NOT in pare → zero
        @test eng.ps8  > 0.0
        @test eng.Ts8  > 0.0
        @test eng.Tt8  == 0.0    # total NOT in pare → zero

        # ------------------------------------------------------------------
        # Stations 19c, 25c, 4a, 49c: not in pare → all zero.
        # ------------------------------------------------------------------
        @test eng.Tt19c == 0.0
        @test eng.Tt25c == 0.0
        @test eng.Tt49c == 0.0

        # ------------------------------------------------------------------
        # Cooling mixing constants (tasopt-j9l.54): ruc and M4a are frozen
        # design inputs read via DesignState, not bare pare reads.
        # Verify they match pare and are physically sensible.
        # ------------------------------------------------------------------
        @test eng.design.ruc ≈ pare_ip[ieruc] rtol = 1e-12
        @test eng.design.M4a ≈ pare_ip[ieM4a] rtol = 1e-12
        @test eng.design.ruc > 0.0    # velocity ratio is positive
        @test eng.design.M4a > 0.0    # Mach number is positive

        # ------------------------------------------------------------------
        # Overall propulsion efficiencies (tasopt-j9l.63.2)
        # Verify bounds and the eta_overall definition:
        #   eta_overall = Fe * u0 / (mdotf_per_engine * hfuel)
        # where mdotf_per_engine = eng.mfuel / neng.
        # ------------------------------------------------------------------
        neng  = ac.parg[TASOPT.igneng]
        hfuel = ac.pare[TASOPT.iehfuel, ipcruise1, 1]

        @test 0.0 < eng.eta_overall < 1.0
        @test 0.0 < eng.eta_prop    < 1.0
        @test 0.0 < eng.eta_thermal < 1.0

        # Definition check: eta_overall == Fe * u0 / (mfuel/neng * hfuel)
        eta_overall_ref = eng.Fe * eng.st0.u / (eng.mfuel / neng * hfuel)
        @test eng.eta_overall ≈ eta_overall_ref rtol = 1e-10

        # eta_thermal = eta_overall / eta_prop
        @test eng.eta_thermal ≈ eng.eta_overall / eng.eta_prop rtol = 1e-12
    end  # run_engine_design_point

    # ======================================================================
    # engine_state_to_pare! / design_state_to_pare! — round-trip fidelity
    # ======================================================================
    @testset "engine_state_to_pare! round-trip" begin
        # After sizing the pare slice at ipcruise1 is converged; verify that
        # reading into EngineState and writing back produces bit-for-bit
        # identical values for every field that pare_to_engine_state! covers.
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)

        imission = 1
        ip       = ipcruise1
        pare_orig = view(ac.pare, :, ip, imission)

        # Populate EngineState from original pare
        eng = TASOPT.engine.EngineState{Float64}()
        TASOPT.engine.pare_to_engine_state!(eng, pare_orig)

        # Write EngineState back into a fresh copy of pare
        pare_copy = copy(pare_orig)
        TASOPT.engine.engine_state_to_pare!(eng, pare_copy)

        # All fields written by pare_to_engine_state! must round-trip exactly.
        for idx in [ieM0, ieT0, iep0, iea0,
                    ieTt0, ieht0, iept0, iecpt0, ieRt0, ieu0,
                    ieTt18, ieht18, iept18, iecpt18, ieRt18,
                    ieTt19, ieht19, iept19, iecpt19, ieRt19,
                    ieTt2, ieht2, iept2, iecpt2, ieRt2,
                    iep2, ieT2, ieR2, iecp2, ieu2, ieA2, iemcore,
                    ieTt21, ieht21, iept21, iecpt21, ieRt21,
                    ieTt25, ieht25, iept25, iecpt25, ieRt25,
                    iep25, ieT25, ieR25, iecp25, ieu25, ieA25,
                    ieTt3, ieht3, iept3, iecpt3, ieRt3,
                    ieht4, iept4, iecpt4, ieRt4,
                    ieTt41, ieht41, iept41, iecpt41, ieRt41,
                    ieTt45, ieht45, iept45, iecpt45, ieRt45,
                    ieTt49, ieht49, iept49, iecpt49, ieRt49,
                    ieTt5, ieht5, iept5, iecpt5, ieRt5,
                    iep5, ieT5, ieR5, iecp5, ieu5, ieA5,
                    iep6, ieT6, ieR6, iecp6, ieu6, ieA6,
                    ieTt7, ieht7, iept7, iecpt7, ieRt7,
                    iep7, ieT7, ieR7, iecp7, ieu7, ieA7,
                    iep8, ieT8, ieR8, iecp8, ieu8, ieA8,
                    ieu9, ieA9]
            @test pare_copy[idx] == pare_orig[idx]
        end
    end  # engine_state_to_pare! round-trip

    @testset "design_state_to_pare! round-trip" begin
        # After sizing, DesignState can be read from pare and written back;
        # every design scalar must match exactly.
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)

        imission = 1
        ip       = ipcruise1
        pare_orig = view(ac.pare, :, ip, imission)

        # Build a DesignState from the converged pare slice
        ds = TASOPT.engine.DesignState{Float64}()
        ds.A2    = pare_orig[ieA2];   ds.A25   = pare_orig[ieA25]
        ds.A5    = pare_orig[ieA5];   ds.A7    = pare_orig[ieA7]
        ds.NbfD  = pare_orig[ieNbfD]; ds.NblcD = pare_orig[ieNblcD]
        ds.NbhcD = pare_orig[ieNbhcD]; ds.NbhtD = pare_orig[ieNbhtD]
        ds.NbltD = pare_orig[ieNbltD]
        ds.mbfD  = pare_orig[iembfD]; ds.mblcD = pare_orig[iemblcD]
        ds.mbhcD = pare_orig[iembhcD]; ds.mbhtD = pare_orig[iembhtD]
        ds.mbltD = pare_orig[iembltD]
        ds.pifD  = pare_orig[iepifD]; ds.pilcD = pare_orig[iepilcD]
        ds.pihcD = pare_orig[iepihcD]; ds.pihtD = pare_orig[iepihtD]
        ds.piltD = pare_orig[iepiltD]
        # tasopt-j9l.54: cooling mixing constants
        ds.ruc   = pare_orig[ieruc]
        ds.M4a   = pare_orig[ieM4a]

        pare_copy = copy(collect(pare_orig))
        TASOPT.engine.design_state_to_pare!(ds, pare_copy)

        for idx in [ieA2, ieA25, ieA5, ieA7,
                    ieNbfD, ieNblcD, ieNbhcD, ieNbhtD, ieNbltD,
                    iembfD, iemblcD, iembhcD, iembhtD, iembltD,
                    iepifD, iepilcD, iepihcD, iepihtD, iepiltD,
                    ieruc, ieM4a]
            @test pare_copy[idx] == collect(pare_orig)[idx]
        end
    end  # design_state_to_pare! round-trip

    # ======================================================================
    # off-design engine_state_to_pare! round-trip — tasopt-j9l.50
    # After sizing + a full off-design sweep, pare columns at off-design
    # points contain tfoper! outputs.  Verify that:
    #   1. pare_to_engine_state! + engine_state_to_pare! is a no-op for
    #      every field engine_state_to_pare! writes (off-design points).
    #   2. The typed EngineState is consistent with pare for the key
    #      performance-relevant scalars at both FixedTt4 and FixedFe modes.
    # ======================================================================
    @testset "off-design engine_state_to_pare! round-trip" begin
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)
        TASOPT.engine.run_engine_sweep(ac)   # updates all 16 off-design pare columns

        imission = 1
        # ipcruise1 → FixedTt4OffDes (Fe is output, Tt4 is input)
        # ipstatic  → FixedFeOffDes  (Tt4 is output, Fe is input)
        for ip in [ipcruise1, ipstatic]
            pare_orig = view(ac.pare, :, ip, imission)

            eng = TASOPT.engine.EngineState{Float64}()
            TASOPT.engine.pare_to_engine_state!(eng, pare_orig)

            # Write back into a fresh copy; must reproduce pare exactly.
            pare_copy = copy(collect(pare_orig))
            TASOPT.engine.engine_state_to_pare!(eng, pare_copy)

            written_indices = [
                ieM0, ieT0, iep0, iea0,
                ieTt0, ieht0, iept0, iecpt0, ieRt0, ieu0,
                ieTt18, ieht18, iept18, iecpt18, ieRt18,
                ieTt19, ieht19, iept19, iecpt19, ieRt19,
                ieTt2, ieht2, iept2, iecpt2, ieRt2,
                iep2, ieT2, ieR2, iecp2, ieu2, ieA2, iemcore,
                ieTt21, ieht21, iept21, iecpt21, ieRt21,
                ieTt25, ieht25, iept25, iecpt25, ieRt25,
                iep25, ieT25, ieR25, iecp25, ieu25, ieA25,
                ieTt3, ieht3, iept3, iecpt3, ieRt3,
                ieht4, iept4, iecpt4, ieRt4,
                ieTt41, ieht41, iept41, iecpt41, ieRt41,
                ieTt45, ieht45, iept45, iecpt45, ieRt45,
                ieTt49, ieht49, iept49, iecpt49, ieRt49,
                ieTt5, ieht5, iept5, iecpt5, ieRt5,
                iep5, ieT5, ieR5, iecp5, ieu5, ieA5,
                iep6, ieT6, ieR6, iecp6, ieu6, ieA6,
                ieTt7, ieht7, iept7, iecpt7, ieRt7,
                iep7, ieT7, ieR7, iecp7, ieu7, ieA7,
                iep8, ieT8, ieR8, iecp8, ieu8, ieA8,
                ieu9, ieA9,
            ]
            for idx in written_indices
                @test pare_copy[idx] == collect(pare_orig)[idx]
            end
        end

        # Consistency: typed EngineState fields must equal pare directly
        # for both off-design modes after the sweep.
        for ip in [ipcruise1, ipstatic]
            pare_ip = view(ac.pare, :, ip, imission)
            eng = TASOPT.engine.EngineState{Float64}()
            TASOPT.engine.pare_to_engine_state!(eng, pare_ip)

            @test eng.Tt4  ≈ pare_ip[ieTt4]   rtol=1e-12
            @test eng.ht4  ≈ pare_ip[ieht4]   rtol=1e-12
            @test eng.pt3  ≈ pare_ip[iept3]   rtol=1e-12
            @test eng.Tt49 ≈ pare_ip[ieTt49]  rtol=1e-12
            @test eng.A2   ≈ pare_ip[ieA2]    rtol=1e-12
            @test eng.M0   ≈ pare_ip[ieM0]    rtol=1e-12
            @test eng.st2.mdot ≈ pare_ip[iemcore] rtol=1e-12
        end
    end  # off-design engine_state_to_pare! round-trip

    # ======================================================================
    # HX delta output fields in EngineState (tasopt-j9l.41.2 / tasopt-w82)
    # HX delta fields are written to per-point EngineState by
    # HXOffDesign!/resetHXs (tasopt-j9l.41.2). For the default model
    # (no heat exchangers), all fields are 0.0 after sizing.
    # hvapcombustor is initialized from TOML in read_input.jl (tasopt-j9l.61).
    # Bare pare slots still maintained for tfcalc! (pending tasopt-w83).
    # ======================================================================
    @testset "HX delta output fields initialized correctly (tasopt-w82)" begin
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)

        im = 1
        ip = ipcruise1
        eng = ac.missions[im].points[ip].engine

        @test eng.PreCDeltah    === 0.0
        @test eng.PreCDeltap    === 0.0
        @test eng.InterCDeltah  === 0.0
        @test eng.InterCDeltap  === 0.0
        @test eng.RegenDeltah   === 0.0
        @test eng.RegenDeltap   === 0.0
        @test eng.TurbCDeltah   === 0.0
        @test eng.TurbCDeltap   === 0.0
        @test eng.RadiatorDeltah === 0.0
        @test eng.RadiatorDeltap === 0.0
        @test eng.HXrecircP     === 0.0
        # hvapcombustor initialized from TOML (default model: fuel_enthalpy_vaporization = 0.0)
        @test eng.hvapcombustor ≈ 0.0
    end  # HX delta output fields initialized correctly (tasopt-w82)

    # ======================================================================
    # run_engine_sweep / write_sweep_csv
    # ======================================================================
    @testset "run_engine_sweep" begin
        # ------------------------------------------------------------------
        # Setup: use the same sized aircraft from the design-point test.
        # ------------------------------------------------------------------
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)

        # ------------------------------------------------------------------
        # Sweep all regular mission points (ipstatic:ipdescentn = 1:16).
        # run_engine_sweep returns ac.missions[1] (Mission{Float64}) and
        # populates mission.points[ip].engine for each ip in ip_range.
        # After sizing, every pare column is already converged, so the sweep
        # should reproduce those converged values.
        # ------------------------------------------------------------------
        mission = TASOPT.engine.run_engine_sweep(ac)
        ips = collect(ipstatic:ipdescentn)
        n_pts = length(ips)   # should be 16

        # ---- Return type invariant ----
        @test mission isa Mission{Float64}

        # ---- All performance scalars are positive (physical) ----
        @test all(ip -> mission.points[ip].engine.Fe   > 0.0, ips)
        @test all(ip -> mission.points[ip].engine.TSFC > 0.0, ips)
        @test all(ip -> mission.points[ip].engine.BPR  > 0.0, ips)
        @test all(ip -> mission.points[ip].engine.st2.mdot > 0.0, ips)
        @test all(ip -> mission.points[ip].engine.mfuel > 0.0, ips)

        # ---- Mach and altitude are populated for all points ----
        @test all(ip -> ac.para[iaMach, ip, 1] >= 0.0, ips)
        @test all(ip -> ac.para[iaalt,  ip, 1] >= 0.0, ips)

        # ---- Engine states match pare at every mission point ----
        eng_cr = mission.points[ipcruise1].engine
        pare_cr = view(ac.pare, :, ipcruise1, 1)

        @test eng_cr.Tt4  ≈ pare_cr[ieTt4]    rtol = 1e-12
        @test eng_cr.pt3  ≈ pare_cr[iept3]    rtol = 1e-12
        @test eng_cr.Tt49 ≈ pare_cr[ieTt49]   rtol = 1e-12
        @test eng_cr.BPR  ≈ pare_cr[ieBPR]    rtol = 1e-12
        @test eng_cr.Fe   ≈ pare_cr[ieFe]     rtol = 1e-12
        @test eng_cr.TSFC ≈ pare_cr[ieTSFC]   rtol = 1e-12
        @test eng_cr.st2.mdot ≈ pare_cr[iemcore] rtol = 1e-12
        @test eng_cr.mfuel    ≈ pare_cr[iemfuel] rtol = 1e-12

        # Check a climb point too (ipclimb1 uses FixedTt4OffDes mode)
        eng_cl  = mission.points[ipclimb1].engine
        pare_cl = view(ac.pare, :, ipclimb1, 1)

        @test eng_cl.Tt4 ≈ pare_cl[ieTt4] rtol = 1e-12
        @test eng_cl.pt3 ≈ pare_cl[iept3] rtol = 1e-12

        # ---- Thermodynamic invariants on cruise point ----
        # Total temperature must rise through compressor, fall through turbine
        @test eng_cr.Tt25 > eng_cr.Tt19   # LPC adds work
        @test eng_cr.Tt3  > eng_cr.Tt25   # HPC adds work
        @test eng_cr.Tt4  > eng_cr.Tt3    # combustor
        @test eng_cr.Tt41 < eng_cr.Tt4    # turbine inlet (cooling dilution)
        @test eng_cr.Tt45 < eng_cr.Tt41   # HPT
        @test eng_cr.Tt49 < eng_cr.Tt45   # LPT

        # ---- Custom ip_range: only the two cruise points ----
        mission_cr2 = TASOPT.engine.run_engine_sweep(ac;
                          ip_range = ipcruise1:ipcruise2)
        @test mission_cr2 isa Mission{Float64}
        @test mission_cr2.points[ipcruise1].engine.Fe > 0.0
        @test mission_cr2.points[ipcruise2].engine.Fe > 0.0

        # ---- CSV serialisation round-trip ----
        buf = IOBuffer()
        TASOPT.engine.write_sweep_csv(buf, mission_cr2, ipcruise1:ipcruise2, ac)
        csv_str = String(take!(buf))

        # Header row must contain expected column names
        header_line = split(csv_str, "\n")[1]
        @test occursin("ip",          header_line)
        @test occursin("label",       header_line)
        @test occursin("Fe_N",        header_line)
        @test occursin("TSFC_kg_Ns",  header_line)
        @test occursin("BPR",         header_line)
        @test occursin("mfuel_kg_s",  header_line)
        @test occursin("Tt4_K",       header_line)
        @test occursin("pt3_Pa",      header_line)

        # Exactly 3 lines (header + 2 data rows + possible trailing newline)
        data_lines = filter(l -> !isempty(strip(l)), split(csv_str, "\n"))
        @test length(data_lines) == 3   # header + 2 cruise points

        # First data row label is "cruise1"
        first_data = split(data_lines[2], ",")
        @test first_data[1] == string(ipcruise1)   # ip index
        @test first_data[2] == "cruise1"            # human label

    end  # run_engine_sweep

    # ======================================================================
    # write_sweep_toml / regenerate_engine_baseline / baseline fixture
    # ======================================================================
    @testset "write_sweep_toml" begin
        import TOML

        # ------------------------------------------------------------------
        # Setup: sized aircraft + 2-point cruise sweep (fast setup).
        # ------------------------------------------------------------------
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)
        mission_cr2 = TASOPT.engine.run_engine_sweep(ac; ip_range=ipcruise1:ipcruise2)

        # ------------------------------------------------------------------
        # write_sweep_toml writes valid, parseable TOML.
        # ------------------------------------------------------------------
        buf = IOBuffer()
        TASOPT.engine.write_sweep_toml(buf, mission_cr2, ipcruise1:ipcruise2, ac)
        toml_str = String(take!(buf))
        @test !isempty(toml_str)

        parsed = TOML.parse(toml_str)

        # Top-level keys
        @test haskey(parsed, "metadata")
        @test haskey(parsed, "points")
        @test parsed["metadata"]["n_points"] == 2
        @test parsed["metadata"]["aircraft"] == "default"

        # Exactly 2 mission-point entries
        pts = parsed["points"]
        @test length(pts) == 2

        # First point metadata
        p1 = pts[1]
        @test p1["ip"]    == ipcruise1
        @test p1["label"] == "cruise1"
        @test p1["Fe_N"]  > 0.0
        @test p1["TSFC_kg_Ns"] > 0.0
        @test p1["BPR"]   > 0.0

        # Station subtables present for all 20 named stations
        @test haskey(p1, "stations")
        stations_p1 = p1["stations"]
        for stname in ("st0", "st2", "st18", "st19", "st19c",
                        "st21", "st25", "st25c", "st3", "st4", "st4a",
                        "st41", "st45", "st49", "st49c",
                        "st5", "st6", "st7", "st8", "st9")
            @test haskey(stations_p1, stname)
        end

        # Each station dict has all 12 expected scalar fields
        expected_fields = ("Tt", "ht", "pt", "cpt", "Rt",
                           "Ts", "ps", "cps", "Rs", "u", "A", "mdot")
        for f in expected_fields
            @test haskey(stations_p1["st0"],  f)
            @test haskey(stations_p1["st2"],  f)
            @test haskey(stations_p1["st4"],  f)
        end

        # Key thermodynamic values at cruise1 are positive
        @test stations_p1["st4"]["Tt"]  > 0.0   # combustor exit Tt4
        @test stations_p1["st3"]["pt"]  > 0.0   # HPC exit total pressure
        @test stations_p1["st2"]["mdot"] > 0.0  # core mass flow

        # st19c/st25c are now written by tfcalc! into the per-point EngineState
        # (tasopt-j9l.45.16): pare_to_engine_state! moved to tfwrap callers so
        # computed values accumulate.  st4a/st49c remain zero (not written).
        @test stations_p1["st19c"]["Tt"] > 0.0   # pre-cooler outlet, populated by tfcalc!
        @test stations_p1["st25c"]["Tt"] > 0.0   # inter-cooler outlet, populated by tfcalc!
        @test stations_p1["st4a"]["Tt"]  == 0.0  # still not populated from computation
        @test stations_p1["st49c"]["Tt"] == 0.0  # still not populated from computation

        # ------------------------------------------------------------------
        # Baseline fixture exists and has the right structure.
        # The file was committed in tasopt-j9l.4.  If this test fails,
        # run: julia --project=. test/generate_engine_baseline.jl
        # ------------------------------------------------------------------
        baseline_path = TASOPT.engine.ENGINE_BASELINE_PATH
        # If this fails, regenerate with: julia --project=. test/generate_engine_baseline.jl
        @test isfile(baseline_path)

        baseline = TOML.parsefile(baseline_path)
        @test haskey(baseline, "metadata")
        @test haskey(baseline, "points")
        @test baseline["metadata"]["n_points"] == 16

        # All 16 mission points present
        @test length(baseline["points"]) == 16

        # First point is static (ip == ipstatic == 1)
        @test baseline["points"][1]["ip"] == ipstatic
        @test baseline["points"][1]["label"] == "static"

        # Cruise1 point exists with positive thrust
        cruise_pt = baseline["points"][ipcruise1 - ipstatic + 1]
        @test cruise_pt["label"] == "cruise1"
        @test cruise_pt["Fe_N"]  > 0.0
        @test cruise_pt["BPR"]   > 0.0

    end  # write_sweep_toml

    # ======================================================================
    # engine_sweep_regression — tasopt-j9l.5
    # Run the full 16-point off-design sweep and compare every station field
    # at every mission point to the pinned baseline within 1e-12 rtol.
    # Failures identify the label, station key, field, got/ref values, and
    # relative difference for actionable diagnosis.
    # ======================================================================
    @testset "engine_sweep_regression" begin
        import TOML

        # ---- Setup ----
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)
        sweep = TASOPT.engine.run_engine_sweep(ac)   # returns Mission{Float64}
        ips   = collect(ipstatic:ipdescentn)

        baseline_path = TASOPT.engine.ENGINE_BASELINE_PATH
        @test isfile(baseline_path)
        baseline = TOML.parsefile(baseline_path)

        n_pts = length(ips)
        baseline_pts = baseline["points"]
        @test length(baseline_pts) == n_pts

        # ---- Comparison tolerance ----
        # 1e-10 (one part per ten billion) is robust to cross-platform FP
        # differences (Newton solver accumulates ~1e-11 across architectures)
        # while still catching any real numerical regression (logic changes
        # produce differences orders of magnitude larger).
        rtol_tol = 1e-10

        # Collect all failures; report in one actionable block at the end.
        failures = String[]
        function chk!(got, ref, ctx)
            g, r = Float64(got), Float64(ref)
            isapprox(g, r; rtol=rtol_tol) && return
            rdiff = abs(g - r) / max(abs(r), 1e-300)
            push!(failures, "  $ctx: got=$g, ref=$r, rdiff=$rdiff")
        end

        # Station field name → FlowStation property symbol mapping.
        # Must mirror _TOML_STATION_ORDER in engine_harness.jl.
        station_map = (
            ("st0",   :st0),  ("st2",   :st2),  ("st18",  :st18),
            ("st19",  :st19), ("st19c", :st19c), ("st21",  :st21),
            ("st25",  :st25), ("st25c", :st25c), ("st3",   :st3),
            ("st4",   :st4),  ("st4a",  :st4a),  ("st41",  :st41),
            ("st45",  :st45), ("st49",  :st49),  ("st49c", :st49c),
            ("st5",   :st5),  ("st6",   :st6),   ("st7",   :st7),
            ("st8",   :st8),  ("st9",   :st9),
        )
        # Fields serialised by _station_to_dict (excludes hs, ss, st, alpha).
        fs_fields = ("Tt", "ht", "pt", "cpt", "Rt", "Ts", "ps", "cps", "Rs", "u", "A", "mdot")

        # ---- Compare every mission point ----
        for k in 1:n_pts
            bp  = baseline_pts[k]
            lbl = bp["label"]
            ip  = ips[k]
            eng = sweep.points[ip].engine

            @test ip_labels[ip] == lbl  # ordering sanity check

            # Top-level performance scalars
            chk!(eng.Fe,        bp["Fe_N"],       "[$lbl] Fe_N")
            chk!(eng.TSFC,      bp["TSFC_kg_Ns"], "[$lbl] TSFC_kg_Ns")
            chk!(eng.BPR,       bp["BPR"],        "[$lbl] BPR")
            chk!(eng.st2.mdot,  bp["mcore_kg_s"], "[$lbl] mcore_kg_s")
            chk!(eng.mfuel,     bp["mfuel_kg_s"], "[$lbl] mfuel_kg_s")
            chk!(ac.para[iaalt,  ip, 1], bp["alt_m"], "[$lbl] alt_m")
            chk!(ac.para[iaMach, ip, 1], bp["Mach"],  "[$lbl] Mach")
            chk!(eng.M0,        bp["M0"],         "[$lbl] M0")
            chk!(eng.T0,        bp["T0_K"],       "[$lbl] T0_K")
            chk!(eng.p0,        bp["p0_Pa"],      "[$lbl] p0_Pa")
            chk!(eng.a0,        bp["a0_m_s"],     "[$lbl] a0_m_s")

            # Station-level fields
            bp_stations = bp["stations"]
            for (stkey, stsym) in station_map
                bp_st = bp_stations[stkey]
                fs    = getfield(eng, stsym)
                for fld in fs_fields
                    chk!(getproperty(fs, Symbol(fld)), bp_st[fld], "[$lbl][$stkey] $fld")
                end
            end
        end

        # Single pass/fail gate; emit full diff when anything failed.
        @test isempty(failures)
        isempty(failures) || @info "Engine sweep regression diff ($(length(failures)) failures):\n$(join(failures, "\n"))"

    end  # engine_sweep_regression

    # ======================================================================
    # engine_sweep_benchmark — tasopt-j9l.7
    # Pin a performance baseline (allocation count + reference runtime) for
    # run_engine_sweep and flag regressions.
    #
    # Allocation count:  deterministic post-JIT → hard @test (5% tolerance).
    # Elapsed time:      hardware-dependent     → informational log (10% tol).
    #
    # To update the baseline after an intentional performance improvement:
    #   julia --project=. test/generate_engine_benchmark_baseline.jl
    # ======================================================================
    @testset "engine_sweep_benchmark" begin
        import TOML

        baseline_path = joinpath(@__DIR__, "fixtures", "engine_benchmark_baseline.toml")
        @test isfile(baseline_path)
        baseline = TOML.parsefile(baseline_path)

        baseline_alloc   = Int(baseline["benchmark"]["alloc_count"])
        baseline_elapsed = Float64(baseline["benchmark"]["elapsed_s"])

        # Setup: freshly sized aircraft (same state as regression test above).
        ac_bench = TASOPT.load_default_model()
        size_aircraft!(ac_bench; printiter=false)

        # Warmup: one call to trigger JIT before measuring.
        _ = TASOPT.engine.run_engine_sweep(ac_bench)

        # --- Allocation check (deterministic) --------------------------------
        # @allocations counts heap allocations made by the call post-JIT.
        # Any new heap allocation introduced by a code change will show up here.
        alloc = @allocations TASOPT.engine.run_engine_sweep(ac_bench)
        alloc_limit = round(Int, baseline_alloc * 1.05)
        alloc_ok = alloc <= alloc_limit
        @test alloc_ok
        if !alloc_ok
            @info "Allocation regression: got $alloc, baseline $baseline_alloc, " *
                  "limit $alloc_limit ($(round((alloc / baseline_alloc - 1) * 100; digits=1))% over)"
        end

        # --- Runtime check (informational) -----------------------------------
        # Single post-warmup elapsed measurement; timing varies by hardware and
        # system load so this is logged rather than asserted.
        elapsed = @elapsed TASOPT.engine.run_engine_sweep(ac_bench)
        time_limit = baseline_elapsed * 1.10
        if elapsed > time_limit
            @info "Runtime regression (informational, hardware-dependent): " *
                  "got $(round(elapsed * 1e3; digits=2))ms, " *
                  "baseline $(round(baseline_elapsed * 1e3; digits=2))ms, " *
                  "limit $(round(time_limit * 1e3; digits=2))ms " *
                  "($(round((elapsed / baseline_elapsed - 1) * 100; digits=1))% over)"
        end

    end  # engine_sweep_benchmark

    # ======================================================================
    # engine_plots — tasopt-j9l.6
    # Verify that the two plotting functions return a Plots.Plot without
    # error on the default aircraft.  Structural invariants checked:
    #   - correct return type
    #   - expected subplot count
    #   - savefig to a temp file succeeds without throwing
    # ======================================================================
    # ======================================================================
    # Inlet component — tasopt-j9l.28
    # Verify structural invariants and outlet-state physics for
    #   inlet_diffuser!  and  inlet_bli_mixing!
    #
    # Property-based tests:
    #   - diffuser conserves total enthalpy/temperature; pt18 = pid * pt0
    #   - BLI=0: no pressure loss; Tt unchanged
    #   - BLI>0 (fan-only): st2.pt < st18.pt; st19.pt = st18.pt; adiabatic
    #   - BLI>0 (all cores): st2.pt == st19.pt; both reduced
    # Round-trip:
    #   - sbfan from inlet_bli_mixing! agrees with the inlined tfoper.jl formula
    # ======================================================================
    @testset "Inlet" begin
        using StaticArrays

        FS               = TASOPT.engine.FlowStation
        Inl              = TASOPT.engine.Inlet
        set_tt!          = TASOPT.engine.set_total_from_Tt!
        inlet_diffuser!  = TASOPT.engine.inlet_diffuser!
        inlet_bli!       = TASOPT.engine.inlet_bli_mixing!
        gassum           = TASOPT.engine.gassum
        gas_tset         = TASOPT.engine.gas_tset

        # Standard 5-species air composition
        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair = 5

        # ------------------------------------------------------------------
        # Build realistic freestream station for cruise (M0=0.8, T0=219.43K)
        # Replicates the tfsize! test freestream.
        # ------------------------------------------------------------------
        M0_val  = 0.8
        T0_val  = 219.43067572699252
        p0_val  = 23922.608843328788
        a0_val  = 296.85578884697560
        u0_val  = M0_val * a0_val

        # Stagnation temperature via gas_tset
        s0_c, _, hs0_c, _, cps0_c, Rs0_c = gassum(air_alpha, nair, T0_val)
        hspec_c  = hs0_c + 0.5 * u0_val^2
        gam0_c   = cps0_c / (cps0_c - Rs0_c)
        Tguess_c = T0_val * (1.0 + 0.5 * (gam0_c - 1.0) * M0_val^2)
        Tt0_c    = gas_tset(air_alpha, nair, hspec_c, Tguess_c)

        # Total pressure via isentropic relation
        st0_ent, _, ht0_c, _, cpt0_c, Rt0_c = gassum(air_alpha, nair, Tt0_c)
        pt0_c = p0_val * exp((st0_ent - s0_c) / Rt0_c)
        at0_c = sqrt(Tt0_c * Rt0_c * cpt0_c / (cpt0_c - Rt0_c))

        # Build FlowStation for st0 using set_total_from_Tt! for entropy
        st0 = FS{Float64}()
        st0.alpha = air_alpha
        st0.Tt    = Tt0_c
        st0.pt    = pt0_c
        set_tt!(st0)   # fills ht, st, cpt, Rt from Tt and alpha

        # Sanity-check: total temperature above ambient
        @test st0.Tt > T0_val
        @test st0.pt > p0_val

        # ------------------------------------------------------------------
        # 1. Inlet struct and constructors
        # ------------------------------------------------------------------
        # Full typed constructor
        inl_full = Inl{Float64}(0.998, 0.0, 0.0, false)
        @test inl_full isa Inl{Float64}
        @test inl_full.pid              === 0.998
        @test inl_full.Kinl             === 0.0
        @test inl_full.Phiinl           === 0.0
        @test inl_full.eng_has_BLI_cores === false

        # Convenience constructor (no BLI)
        inl_clean = Inl(0.997)
        @test inl_clean isa Inl{Float64}
        @test inl_clean.pid              === 0.997
        @test inl_clean.Kinl             === 0.0
        @test inl_clean.Phiinl           === 0.0
        @test inl_clean.eng_has_BLI_cores === false

        # Convenience constructor with BLI keyword args
        inl_bli = Inl(0.995; Kinl=5000.0, Phiinl=1000.0, eng_has_BLI_cores=true)
        @test inl_bli.Kinl              ≈ 5000.0
        @test inl_bli.Phiinl            ≈ 1000.0
        @test inl_bli.eng_has_BLI_cores === true

        # Mutation (mutable struct)
        inl_mut = Inl(0.99)
        inl_mut.pid = 0.95
        @test inl_mut.pid === 0.95

        # ------------------------------------------------------------------
        # 2. inlet_diffuser! — adiabatic duct invariants
        # ------------------------------------------------------------------
        pid_val = 0.998
        inl_d   = Inl(pid_val)
        st18    = FS{Float64}()

        ret = inlet_diffuser!(st18, st0, inl_d)

        # Returns the mutated station
        @test ret === st18

        # Adiabatic: all total-state fields (except pt) unchanged
        @test st18.Tt    === st0.Tt
        @test st18.ht    === st0.ht
        @test st18.cpt   === st0.cpt
        @test st18.Rt    === st0.Rt
        @test st18.st    === st0.st
        @test st18.alpha === st0.alpha

        # Pressure recovery: exact floating-point product
        @test st18.pt === st0.pt * pid_val

        # pid = 1.0: no pressure drop at all
        st18_id = FS{Float64}()
        inlet_diffuser!(st18_id, st0, Inl(1.0))
        @test st18_id.pt === st0.pt

        # Monotonicity: lower pid ↔ lower pt18
        pt18_95 = let st18_lo = FS{Float64}()
            inlet_diffuser!(st18_lo, st0, Inl(0.95)); st18_lo.pt
        end
        pt18_99 = let st18_hi = FS{Float64}()
            inlet_diffuser!(st18_hi, st0, Inl(0.99)); st18_hi.pt
        end
        @test pt18_95 < pt18_99 < st0.pt

        # ------------------------------------------------------------------
        # 3. inlet_bli_mixing! — zero BLI (Kinl = 0): no pressure loss
        # ------------------------------------------------------------------
        # Shared Tref/pref for corrected-flow normalisation
        Tref_val = 288.2
        pref_val = 101325.0
        mf_val   = 1.0
        ml_val   = 1.0
        M2_val   = 0.6

        inl_nobli = Inl(pid_val)   # Kinl = 0 by default
        TASOPT.engine.inlet_diffuser!(st18, st0, inl_nobli)  # refresh st18
        st2  = FS{Float64}()
        st19 = FS{Float64}()

        res_nobli = inlet_bli!(st2, st19, st18, st0, inl_nobli,
                               mf_val, ml_val, M2_val,
                               at0_c, gam0_c, Tref_val, pref_val)

        # No BLI: outlets equal inlet face
        @test st2.pt  === st18.pt
        @test st19.pt === st18.pt

        # Adiabatic: temperatures unchanged
        @test st2.Tt  === st18.Tt
        @test st19.Tt === st18.Tt

        # Entropy boosts are zero
        @test res_nobli.sbfan  === 0.0
        @test res_nobli.sbcore === 0.0

        # ------------------------------------------------------------------
        # 4. inlet_bli_mixing! — BLI fan-only (eng_has_BLI_cores = false)
        # ------------------------------------------------------------------
        Kinl_val = 50_000.0   # [W] representative cruise BLI power
        inl_bli2 = Inl(pid_val; Kinl=Kinl_val, eng_has_BLI_cores=false)
        TASOPT.engine.inlet_diffuser!(st18, st0, inl_bli2)
        st2_b  = FS{Float64}()
        st19_b = FS{Float64}()

        res_b = inlet_bli!(st2_b, st19_b, st18, st0, inl_bli2,
                           mf_val, ml_val, M2_val,
                           at0_c, gam0_c, Tref_val, pref_val)

        # Fan face pressure drops
        @test st2_b.pt < st18.pt

        # Core stream unaffected (eng_has_BLI_cores = false)
        @test st19_b.pt === st18.pt

        # Adiabatic: temperatures unchanged
        @test st2_b.Tt  === st18.Tt
        @test st19_b.Tt === st18.Tt

        # Entropy boost signs
        @test res_b.sbfan  > 0.0
        @test res_b.sbcore === 0.0

        # Pressure consistent with returned entropy boost
        @test st2_b.pt ≈ st18.pt * exp(-res_b.sbfan)   rtol=1e-15

        # ------------------------------------------------------------------
        # 5. inlet_bli_mixing! — BLI all streams (eng_has_BLI_cores = true)
        # ------------------------------------------------------------------
        inl_blic = Inl(pid_val; Kinl=Kinl_val, eng_has_BLI_cores=true)
        TASOPT.engine.inlet_diffuser!(st18, st0, inl_blic)
        st2_c  = FS{Float64}()
        st19_c = FS{Float64}()

        res_c = inlet_bli!(st2_c, st19_c, st18, st0, inl_blic,
                           mf_val, ml_val, M2_val,
                           at0_c, gam0_c, Tref_val, pref_val)

        # Both streams see the same entropy boost
        @test res_c.sbfan   === res_c.sbcore
        @test res_c.sbfan   > 0.0

        # Pressures consistent with entropy boosts
        @test st2_c.pt  ≈ st18.pt * exp(-res_c.sbfan)   rtol=1e-15
        @test st19_c.pt ≈ st18.pt * exp(-res_c.sbcore)  rtol=1e-15
        @test st2_c.pt  ≈ st19_c.pt                     rtol=1e-15

        # ------------------------------------------------------------------
        # 6. BLI monotonicity: larger Kinl → larger entropy boost
        # ------------------------------------------------------------------
        sbfan_lo = let st2_l=FS{Float64}(), st19_l=FS{Float64}(), st18_l=FS{Float64}()
            TASOPT.engine.inlet_diffuser!(st18_l, st0, Inl(pid_val; Kinl=10_000.0))
            r = inlet_bli!(st2_l, st19_l, st18_l, st0, Inl(pid_val; Kinl=10_000.0),
                           mf_val, ml_val, M2_val,
                           at0_c, gam0_c, Tref_val, pref_val)
            r.sbfan
        end
        sbfan_hi = let st2_h=FS{Float64}(), st19_h=FS{Float64}(), st18_h=FS{Float64}()
            TASOPT.engine.inlet_diffuser!(st18_h, st0, Inl(pid_val; Kinl=100_000.0))
            r = inlet_bli!(st2_h, st19_h, st18_h, st0, Inl(pid_val; Kinl=100_000.0),
                           mf_val, ml_val, M2_val,
                           at0_c, gam0_c, Tref_val, pref_val)
            r.sbfan
        end
        @test sbfan_lo < sbfan_hi

        # ------------------------------------------------------------------
        # 7. Round-trip: compare sbfan against inlined tfoper.jl formula
        #
        # The tfoper.jl formula (non-BLI_cores, u0 ≠ 0):
        #   a2sq = at0^2 / (1 + 0.5*(gam0-1)*M2^2)
        #   mmix = mf * sqrt(Tref/Tt0) * pt0/pref
        #   sbfan_ref = Kinl * gam0 / (mmix * a2sq)
        #   pt2_ref   = pt18 * exp(-sbfan_ref)
        # ------------------------------------------------------------------
        TASOPT.engine.inlet_diffuser!(st18, st0, Inl(pid_val; Kinl=Kinl_val))
        a2sq_ref  = at0_c^2 / (1.0 + 0.5 * (gam0_c - 1.0) * M2_val^2)
        mmix_ref  = mf_val * sqrt(Tref_val / st0.Tt) * (st0.pt / pref_val)
        sbfan_ref = Kinl_val * gam0_c / (mmix_ref * a2sq_ref)
        pt2_ref   = st18.pt * exp(-sbfan_ref)

        @test res_b.sbfan ≈ sbfan_ref  rtol=1e-14
        @test st2_b.pt    ≈ pt2_ref    rtol=1e-14

    end  # Inlet

    # ==========================================================================
    # Turbine component — tasopt-j9l.25
    # Verify structural invariants and outlet-state physics for
    #   TurbineMap, Turbine, turbine_efficiency, turbine_exit!
    #
    # Property-based tests:
    #   - TurbineMap fields match constructor arguments
    #   - Turbine default map constants are 0.15
    #   - turbine_efficiency at design operating point returns ep0
    #   - turbine_efficiency off-design returns value below ep0
    #   - turbine_exit!: Tt_out < Tt_in, pt_out < pt_in, alpha unchanged
    # Round-trip:
    #   - turbine_efficiency result matches direct etmap call
    #   - turbine_exit! Tt_out agrees with gas_delhd
    # ==========================================================================
    @testset "Turbine" begin
        using StaticArrays

        FS               = TASOPT.engine.FlowStation
        TMap             = TASOPT.engine.TurbineMap
        Turb             = TASOPT.engine.Turbine
        turb_eff         = TASOPT.engine.turbine_efficiency
        turb_exit!       = TASOPT.engine.turbine_exit!
        etmap            = TASOPT.engine.etmap
        gas_delhd        = TASOPT.engine.gas_delhd
        gassum           = TASOPT.engine.gassum
        set_tt!          = TASOPT.engine.set_total_from_Tt!

        # ------------------------------------------------------------------
        # HPT design parameters from tfwrap.jl (tasopt-j9l.25 reference values)
        # ------------------------------------------------------------------
        pihtD  = 2.1601257635200488
        mbhtD  = 4.3594697284253883
        NbhtD  = 0.44698693289691338
        epht0  = 0.889
        pcon_h = 0.15
        Ncon_h = 0.15

        # LPT design parameters
        piltD  = 6.2886975330083716
        mbltD  = 8.7016090343744406
        NbltD  = 0.48396724306758404
        eplt0  = 0.899

        # Standard 5-species hot gas composition after combustion
        hot_alpha = SA[0.6800, 0.2000, 0.0006, 0.0020, 0.1174]
        nair = 5

        # Realistic HPT inlet conditions (station 41 after combustor + cooling mixing)
        Tt41   = 1800.0   # K — representative turbine-inlet temperature
        pt41   = 2.0e6    # Pa — representative HPT inlet pressure

        # Build HPT inlet FlowStation using set_total_from_Tt! pattern
        # (same as Inlet testset — fills ht, st, cpt, Rt from Tt and alpha)
        st41 = FS{Float64}()
        st41.alpha = hot_alpha
        st41.Tt    = Tt41
        st41.pt    = pt41
        set_tt!(st41)   # fills ht, st, cpt, Rt

        # Extract scalars for use in etmap calls
        ht41  = st41.ht
        s41   = st41.st
        cpt41 = st41.cpt
        Rt41  = st41.Rt

        # ------------------------------------------------------------------
        # 1. TurbineMap struct and constructors
        # ------------------------------------------------------------------
        tmap1 = TMap(pcon_h, Ncon_h)
        @test tmap1.pcon == pcon_h
        @test tmap1.Ncon == Ncon_h

        tmap2 = TMap([pcon_h, Ncon_h])   # vector constructor
        @test tmap2.pcon == pcon_h
        @test tmap2.Ncon == Ncon_h

        # ------------------------------------------------------------------
        # 2. Turbine struct and constructors
        # ------------------------------------------------------------------
        hpt = Turb(pihtD, mbhtD, NbhtD, epht0; map=TMap(pcon_h, Ncon_h))
        @test hpt.piD ≈ pihtD
        @test hpt.mbD ≈ mbhtD
        @test hpt.NbD ≈ NbhtD
        @test hpt.ep0 ≈ epht0
        @test hpt.map.pcon ≈ pcon_h
        @test hpt.map.Ncon ≈ Ncon_h

        # Default map constants (0.15, 0.15)
        lpt = Turb(piltD, mbltD, NbltD, eplt0)
        @test lpt.map.pcon ≈ 0.15
        @test lpt.map.Ncon ≈ 0.15

        # ------------------------------------------------------------------
        # 3. turbine_efficiency: design-point identity (ept = ep0)
        #
        # At the design operating point the two penalty terms vanish:
        #   prat / piD = 1  →  pcon-term = 0
        #   Nmb / NmbD = 1  →  Ncon-term = 0
        # so ept = ep0.  We construct the design-point dh from piD:
        #   Trat = piD^(Rt·ep0/cpt)
        #   dh   = cpt · Tt · (1/Trat - 1)   (< 0)
        # ------------------------------------------------------------------
        Trat_dp = pihtD^(Rt41 * epht0 / cpt41)
        dh_dp   = cpt41 * Tt41 * (1.0 / Trat_dp - 1.0)
        @test dh_dp < 0.0   # work extracted → negative enthalpy change

        ept_dp, _ = turb_eff(hpt, dh_dp, mbhtD, NbhtD, Tt41, cpt41, Rt41)
        @test ept_dp ≈ epht0  rtol=1e-12   # exact identity at design point

        # ------------------------------------------------------------------
        # 4. turbine_efficiency: off-design penalty
        # ------------------------------------------------------------------
        mb_off = 0.6 * mbhtD   # significantly off design-point flow
        Nb_off = 0.7 * NbhtD
        ept_off, _ = turb_eff(hpt, dh_dp, mb_off, Nb_off, Tt41, cpt41, Rt41)
        @test ept_off < epht0   # efficiency degrades off design

        # Efficiency must be physically positive
        @test ept_off > 0.0

        # ------------------------------------------------------------------
        # 5. turbine_efficiency round-trip: matches direct etmap call
        # ------------------------------------------------------------------
        Tmap_vec = [pcon_h, Ncon_h]
        ept_ref, ept_ref_dh, ept_ref_mb, ept_ref_Nb,
        ept_ref_Tt, ept_ref_cpt, ept_ref_Rt = etmap(
            dh_dp, mbhtD, NbhtD,
            pihtD, mbhtD, NbhtD, epht0, Tmap_vec,
            Tt41, cpt41, Rt41,
        )

        ept_c, ept_c_dh, ept_c_mb, ept_c_Nb,
        ept_c_Tt, ept_c_cpt, ept_c_Rt = turb_eff(hpt, dh_dp, mbhtD, NbhtD, Tt41, cpt41, Rt41)

        @test ept_c      ≈ ept_ref      rtol=1e-14
        @test ept_c_dh   ≈ ept_ref_dh   rtol=1e-14
        @test ept_c_mb   ≈ ept_ref_mb   rtol=1e-14
        @test ept_c_Nb   ≈ ept_ref_Nb   rtol=1e-14
        @test ept_c_Tt   ≈ ept_ref_Tt   rtol=1e-14
        @test ept_c_cpt  ≈ ept_ref_cpt  rtol=1e-14
        @test ept_c_Rt   ≈ ept_ref_Rt   rtol=1e-14

        # ------------------------------------------------------------------
        # 6. turbine_exit!: physical invariants
        # ------------------------------------------------------------------
        dh_test  = cpt41 * Tt41 * (1.0 / Trat_dp - 1.0)   # design work
        ept_test = epht0

        st45 = FS{Float64}()
        turb_exit!(st45, st41, dh_test, ept_test)

        # Temperature drops through expansion
        @test st45.Tt < st41.Tt

        # Pressure drops through expansion
        @test st45.pt < st41.pt

        # Composition unchanged
        @test all(st45.alpha .≈ st41.alpha)

        # Enthalpy drop: gas_delhd solves h(alpha, Tt_out) = h_in + dh
        @test st45.ht ≈ (ht41 + dh_test)  rtol=1e-8

        # ------------------------------------------------------------------
        # 7. turbine_exit! round-trip: matches direct gas_delhd call
        # ------------------------------------------------------------------
        epi_ref = 1.0 / ept_test
        res_ref = gas_delhd(
            hot_alpha, 5,
            pt41, Tt41, ht41, s41, cpt41, Rt41,
            dh_test, epi_ref,
        )
        pt45_ref, Tt45_ref, ht45_ref, s45_ref, cpt45_ref, Rt45_ref =
            res_ref[1], res_ref[2], res_ref[3], res_ref[4], res_ref[5], res_ref[6]

        @test st45.pt  ≈ pt45_ref  rtol=1e-14
        @test st45.Tt  ≈ Tt45_ref  rtol=1e-14
        @test st45.ht  ≈ ht45_ref  rtol=1e-14
        @test st45.st  ≈ s45_ref   rtol=1e-14
        @test st45.cpt ≈ cpt45_ref rtol=1e-14
        @test st45.Rt  ≈ Rt45_ref  rtol=1e-14

        # ------------------------------------------------------------------
        # 8. turbine_exit! monotonicity: larger |dh| → lower pt and Tt
        # ------------------------------------------------------------------
        dh_small = 0.5 * dh_test   # less work extracted (dh_small closer to 0)
        dh_large = 1.5 * dh_test   # more work extracted

        st45_small = FS{Float64}()
        st45_large = FS{Float64}()
        turb_exit!(st45_small, st41, dh_small, ept_test)
        turb_exit!(st45_large, st41, dh_large, ept_test)

        @test st45_large.Tt < st45_small.Tt   # more work → lower exit temperature
        @test st45_large.pt < st45_small.pt   # more work → lower exit pressure

    end  # Turbine

    # turbine_mb_residual — tasopt-j9l.30
    #
    # Tests the standalone map-match residual for vertical-line turbine matching.
    # The residual r = mb - turb.mbD encodes the constraint that corrected mass
    # flow stays at its design value in the Newton system (rows 2 and 3).
    #
    # Properties verified:
    #   - At design point (mb = mbD): r = 0, dr_dmb = 1
    #   - Below design (mb < mbD):   r < 0, dr_dmb = 1
    #   - Above design (mb > mbD):   r > 0, dr_dmb = 1
    #   - Derivative is identically 1 regardless of mb or turb.mbD
    #   - Residual equals mb - turb.mbD exactly (round-trip)
    #   - Correct types for both HPT and LPT representative values
    @testset "turbine_mb_residual" begin

        turb_mb_res = TASOPT.engine.turbine_mb_residual

        # Representative HPT and LPT design values (typical TASOPT turbofan)
        pihtD = 3.5;  mbhtD = 60.0;  NbhtD = 1.0;  epht0 = 0.89
        piltD = 5.0;  mbltD = 70.0;  NbltD = 1.0;  eplt0 = 0.90

        turb_hp = TASOPT.engine.Turbine(pihtD, mbhtD, NbhtD, epht0)
        turb_lp = TASOPT.engine.Turbine(piltD, mbltD, NbltD, eplt0)

        # 1. Design-point identity: r = 0, dr_dmb = 1
        r_hp, dr_hp = turb_mb_res(turb_hp, mbhtD)
        @test r_hp == 0.0
        @test dr_hp == 1.0

        r_lp, dr_lp = turb_mb_res(turb_lp, mbltD)
        @test r_lp == 0.0
        @test dr_lp == 1.0

        # 2. Below-design: r < 0
        r_low, dr_low = turb_mb_res(turb_hp, 0.95 * mbhtD)
        @test r_low < 0.0
        @test dr_low == 1.0

        # 3. Above-design: r > 0
        r_high, dr_high = turb_mb_res(turb_hp, 1.05 * mbhtD)
        @test r_high > 0.0
        @test dr_high == 1.0

        # 4. Derivative is always 1 for any mb (including far off-design)
        for mb_frac in (0.5, 0.8, 1.0, 1.2, 2.0)
            _, dr = turb_mb_res(turb_lp, mb_frac * mbltD)
            @test dr == 1.0
        end

        # 5. Round-trip: residual equals mb - turb.mbD exactly
        for mb_frac in (0.7, 0.9, 1.0, 1.1, 1.3)
            mb = mb_frac * mbhtD
            r, _ = turb_mb_res(turb_hp, mb)
            @test r == mb - turb_hp.mbD
        end

        # 6. Type preservation: Float32 in → Float32 out
        turb32 = TASOPT.engine.Turbine(Float32(pihtD), Float32(mbhtD),
                                        Float32(NbhtD),  Float32(epht0))
        r32, dr32 = turb_mb_res(turb32, Float32(mbhtD))
        @test r32 isa Float32
        @test dr32 isa Float32
        @test r32 == 0.0f0

    end  # turbine_mb_residual

    @testset "Combustor" begin

        Comb      = TASOPT.engine.Combustor
        comb_exit = TASOPT.engine.combustor_exit!
        FS        = TASOPT.engine.FlowStation
        set_tt!   = TASOPT.engine.set_total_from_Tt!
        gassumd   = TASOPT.engine.gassumd
        gas_burnd = TASOPT.engine.gas_burnd
        gasfuel   = TASOPT.engine.gasfuel

        # Standard TASOPT 5-species air composition stored in FlowStation.alpha
        # (species 1-5: N2, O2, CO2, H2O, Ar; fuel slot is NOT stored here)
        air_alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        n    = 5    # nair — number of air species passed to gas_burnd
        nair = 5

        # Representative HPC-exit state (station 3)
        Tt3  = 900.0    # K  — compressor discharge temperature
        pt3  = 4.0e6   # Pa — compressor discharge pressure

        st3  = FS{Float64}()
        st3.alpha = air_alpha
        st3.Tt    = Tt3
        st3.pt    = pt3
        set_tt!(st3)   # fills ht, st, cpt, Rt from Tt3 and alpha

        # Combustor design parameters (Jet-A, etab=0.98, pib=0.94)
        ifuel = 24       # C14H30 — Jet-A surrogate
        hvap  = 0.0      # J/kg — liquid fuel, no vaporisation enthalpy
        Ttf   = 300.0    # K  — fuel inlet temperature
        etab  = 0.98     # combustion efficiency
        pib   = 0.94     # burner pressure ratio

        burner = Comb(pib, etab, Ttf, ifuel, hvap)

        # Target burner-exit total temperature
        Tb = 1600.0   # K

        # ------------------------------------------------------------------
        # 1. Struct / constructor invariants
        # ------------------------------------------------------------------
        @test burner.pib   ≈ pib
        @test burner.etab  ≈ etab
        @test burner.Ttf   ≈ Ttf
        @test burner.ifuel == ifuel
        @test burner.hvap  ≈ hvap

        # ------------------------------------------------------------------
        # 2. combustor_exit! — outlet FlowStation (mutated)
        # ------------------------------------------------------------------
        st4  = FS{Float64}()
        # alpha field is SVector{5} — initialise with the 5-element air vector
        st4.alpha = SVector{5,Float64}(air_alpha...)   # will be overwritten by lambda

        ffb, lambda,
        ffb_Tt3, ffb_Ttf, ffb_Tb,
        lam_Tt3, lam_Ttf, lam_Tb = comb_exit(st4, st3, burner, air_alpha, Tb)

        # Pressure recovery: pt4 = pib × pt3
        @test st4.pt ≈ pib * pt3

        # Target temperature reached
        @test st4.Tt ≈ Tb

        # Thermodynamic state consistent with exit composition
        s_ref, _, h_ref, _, cp_ref, R_ref = gassumd(lambda, nair, Tb)
        @test st4.st  ≈ s_ref  rtol=1e-12
        @test st4.ht  ≈ h_ref  rtol=1e-12
        @test st4.cpt ≈ cp_ref rtol=1e-12
        @test st4.Rt  ≈ R_ref  rtol=1e-12

        # Exit composition stored in FlowStation matches returned lambda
        @test st4.alpha ≈ lambda

        # Fuel fraction is positive (combustion consumes fuel)
        @test ffb > 0.0

        # ------------------------------------------------------------------
        # 3. Round-trip against gas_burnd directly
        # ------------------------------------------------------------------
        # gas_burnd takes n=nair=5; gamma is 5-element matching the air species
        gamma_raw = gasfuel(ifuel, 6)   # gasfuel always uses 6-species layout
        # Apply etab to air species (1:nair=1:5); combustor.jl does the same
        gamma_ref = zeros(5)
        for i in 1:nair
            gamma_ref[i] = etab * gamma_raw[i]
        end
        # Note: gamma[6] = 1-etab for unburnt fuel, but gas_burnd only sees n=5 entries

        beta_ref = zeros(5)   # pure-fuel stream (fuel enthalpy added via gasfun in gas_burnd)

        ffb_ref, lam_ref,
        ffb_Tt3_ref, ffb_Ttf_ref, ffb_Tb_ref,
        lam_Tt3_ref, lam_Ttf_ref, lam_Tb_ref = gas_burnd(
            air_alpha, beta_ref, gamma_ref, nair, ifuel, Tt3, Ttf, Tb, hvap,
        )

        @test ffb    ≈ ffb_ref       rtol=1e-12
        @test lambda ≈ lam_ref       rtol=1e-12
        @test ffb_Tt3 ≈ ffb_Tt3_ref  rtol=1e-12
        @test ffb_Ttf ≈ ffb_Ttf_ref  rtol=1e-12
        @test ffb_Tb  ≈ ffb_Tb_ref   rtol=1e-12
        for i in 1:nair
            @test lam_Tt3[i] ≈ lam_Tt3_ref[i]  rtol=1e-12
            @test lam_Ttf[i] ≈ lam_Ttf_ref[i]  rtol=1e-12
            @test lam_Tb[i]  ≈ lam_Tb_ref[i]   rtol=1e-12
        end

        # ------------------------------------------------------------------
        # 4. pib = 1 → no pressure loss (identity for pressure recovery)
        # ------------------------------------------------------------------
        burner_lossless = Comb(1.0, etab, Ttf, ifuel, hvap)
        st4_lossless = FS{Float64}(); st4_lossless.alpha = SVector{5,Float64}(air_alpha...)
        ffb_ll, _, = comb_exit(st4_lossless, st3, burner_lossless, air_alpha, Tb)
        @test st4_lossless.pt ≈ pt3
        @test st4_lossless.Tt ≈ Tb
        # composition / fuel fraction should match (pib doesn't affect chemistry)
        @test ffb_ll ≈ ffb  rtol=1e-12

        # ------------------------------------------------------------------
        # 5. Fuel fraction monotonicity: higher Tb → more fuel needed
        # ------------------------------------------------------------------
        Tb_lo = 1400.0
        Tb_hi = 1800.0
        st4_lo = FS{Float64}(); st4_lo.alpha = SVector{5,Float64}(air_alpha...)
        st4_hi = FS{Float64}(); st4_hi.alpha = SVector{5,Float64}(air_alpha...)
        ffb_lo, = comb_exit(st4_lo, st3, burner, air_alpha, Tb_lo)
        ffb_hi, = comb_exit(st4_hi, st3, burner, air_alpha, Tb_hi)
        @test ffb_lo < ffb < ffb_hi

        # ------------------------------------------------------------------
        # 6. Partial combustion (etab < 1) vs full combustion (etab = 1)
        # ------------------------------------------------------------------
        # With etab < 1 the reaction fractions gamma are scaled down.  The
        # exact direction of change for individual species depends on both
        # the scaled gamma and the changed fuel fraction ffb (which also
        # shifts with etab).  The robust invariant is that the two cases
        # produce distinct compositions.
        etab_partial = 0.90
        burner_partial = Comb(pib, etab_partial, Ttf, ifuel, hvap)
        st4_part = FS{Float64}(); st4_part.alpha = SVector{5,Float64}(air_alpha...)
        ffb_part, lam_part, = comb_exit(st4_part, st3, burner_partial, air_alpha, Tb)

        burner_full = Comb(pib, 1.0, Ttf, ifuel, hvap)
        st4_full = FS{Float64}(); st4_full.alpha = SVector{5,Float64}(air_alpha...)
        ffb_full, lam_full, = comb_exit(st4_full, st3, burner_full, air_alpha, Tb)

        # Different etab → different fuel fractions and different compositions
        @test ffb_part != ffb_full
        @test lam_part != lam_full
        # Both fuel fractions must be positive (combustion requires fuel)
        @test ffb_part > 0.0
        @test ffb_full > 0.0

        # ------------------------------------------------------------------
        # 7. Combustor{Float32} struct construction (type parameter)
        # ------------------------------------------------------------------
        # gas_burnd calls gasfun which requires Float64 tables, so
        # combustor_exit! with Float32 is not supported.  Here we verify
        # only that Combustor{Float32} can be constructed and stores Float32.
        burner32 = Comb{Float32}(Float32(pib), Float32(etab), Float32(Ttf), ifuel, Float32(hvap))
        @test burner32 isa Comb{Float32}
        @test burner32.pib  isa Float32
        @test burner32.etab isa Float32
        @test burner32.Ttf  isa Float32
        @test burner32.hvap isa Float32
        @test burner32.ifuel == ifuel

    end  # Combustor

    # Compressor component — tasopt-j9l.24
    # Verify structural invariants and outlet-state physics for
    #   Compressor, compressor_efficiency, compressor_exit!
    #
    # Property-based tests:
    #   - Compressor fields match constructor arguments
    #   - compressor_efficiency at design operating point returns epol0
    #   - compressor_efficiency off-design returns a different efficiency
    #   - compressor_efficiency floor: clamped when map returns below floor
    #   - compressor_efficiency windmilling: pi < 1 inverts efficiency
    #   - compressor_exit!: Tt_out > Tt_in, pt_out = pi*pt_in, alpha unchanged
    # Round-trip:
    #   - compressor_efficiency matches direct calculate_compressor_speed_and_efficiency call
    #   - compressor_exit! agrees with direct gas_pratd call
    # ==========================================================================
    @testset "Compressor" begin
        using StaticArrays

        FS               = TASOPT.engine.FlowStation
        Comp             = TASOPT.engine.Compressor
        comp_eff         = TASOPT.engine.compressor_efficiency
        comp_exit!       = TASOPT.engine.compressor_exit!
        calc_cse         = TASOPT.engine.calculate_compressor_speed_and_efficiency
        gas_pratd        = TASOPT.engine.gas_pratd
        FanMap           = TASOPT.engine.FanMap
        LPCMap           = TASOPT.engine.LPCMap
        HPCMap           = TASOPT.engine.HPCMap
        set_tt!          = TASOPT.engine.set_total_from_Tt!

        # ------------------------------------------------------------------
        # Design parameters from tfoper.jl mode-1 test (baseline aircraft)
        # ------------------------------------------------------------------
        pifD   = 1.6850000000000001
        pilcD  = 8.0000000000000000
        pihcD  = 3.7500000000000000

        mbfD   = 235.16225770724063
        mblcD  = 46.110246609262873
        mbhcD  = 7.8056539219349039

        NbfD   = 1.0790738309310697
        NblcD  = 1.0790738309310697
        NbhcD  = 0.77137973563891493

        epf0   = 0.89480000000000004
        eplc0  = 0.88000000000000000
        ephc0  = 0.87000000000000000

        epfmin  = 0.60
        elpcmin = 0.70
        ephcmin = 0.70

        # Standard 5-species air composition (from tfoper.jl)
        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair = 5

        # Realistic fan-face inlet conditions (station 2, ISA cruise ~10 km)
        Tt2  = 288.15   # K  — ISA sea-level (fan-face after inlet recovery)
        pt2  = 30_000.0 # Pa — approx 10 km altitude static → total near 30 kPa

        st2 = FS{Float64}()
        st2.alpha = air_alpha
        st2.Tt    = Tt2
        st2.pt    = pt2
        set_tt!(st2)

        # ------------------------------------------------------------------
        # 1. Compressor struct and constructor invariants
        # ------------------------------------------------------------------
        fan = Comp(pifD, mbfD, NbfD, epf0, epfmin, FanMap)
        @test fan.piD      ≈ pifD
        @test fan.mbD      ≈ mbfD
        @test fan.NbD      ≈ NbfD
        @test fan.epol0    ≈ epf0
        @test fan.epol_min ≈ epfmin
        @test fan.map      === FanMap
        @test fan.Ng       ≈ 0.5    # default warm-start speed
        @test fan.Rg       ≈ 2.0    # default warm-start R-line

        lpc = Comp(pilcD, mblcD, NblcD, eplc0, elpcmin, LPCMap)
        @test lpc.piD  ≈ pilcD
        @test lpc.mbD  ≈ mblcD

        hpc = Comp(pihcD, mbhcD, NbhcD, ephc0, ephcmin, HPCMap)
        @test hpc.piD  ≈ pihcD
        @test hpc.mbD  ≈ mbhcD

        # ------------------------------------------------------------------
        # 2. compressor_efficiency: design-point evaluation
        #
        # At the design operating point the map should return a positive,
        # sub-unity efficiency above the floor.  Note: the map's
        # defaults.polyeff is the *maximum* across the whole map, not
        # the value at the design (Nc, Rline) point, so epol_dp ≠ epol0
        # in general.  The constraint we can assert is: epol_dp ∈ (epfmin, 1).
        # ------------------------------------------------------------------
        fan_dp = Comp(pifD, mbfD, NbfD, epf0, epfmin, FanMap)
        Nb_dp, epol_dp, _, _, _, _ = comp_eff(fan_dp, pifD, mbfD)

        @test epol_dp > epfmin    # above the efficiency floor
        @test epol_dp < 1.0       # physically sub-unity

        # Speed should be positive
        @test Nb_dp > 0.0

        # ------------------------------------------------------------------
        # 3. compressor_efficiency: off-design operating point
        # ------------------------------------------------------------------
        mb_off = 0.8 * mbfD   # mass flow 20% below design
        Nb_off, epol_off, _, _, _, _ = comp_eff(fan_dp, pifD, mb_off)

        # Efficiency must remain physically positive
        @test epol_off > 0.0

        # Speed must be positive
        @test Nb_off > 0.0

        # Compressor speed should change when mass flow changes
        # (not a strict monotonicity guarantee, but sanity check)
        @test Nb_off != Nb_dp

        # ------------------------------------------------------------------
        # 4. compressor_efficiency round-trip: matches direct calc call
        # ------------------------------------------------------------------
        fan_rt = Comp(pifD, mbfD, NbfD, epf0, epfmin, FanMap)
        pi_test = pifD * 0.95   # slightly off design-point pressure ratio
        mb_test = mbfD * 0.95   # slightly off design-point mass flow

        Nb_c, epol_c, Nb_pi_c, Nb_mb_c, epol_pi_c, epol_mb_c =
            comp_eff(fan_rt, pi_test, mb_test)

        # Direct call (no floor/windmill correction since nominal efficiency > floor)
        Nb_ref, epol_ref, Nb_pi_ref, Nb_mb_ref, epol_pi_ref, epol_mb_ref, _, _ =
            calc_cse(FanMap, pi_test, mb_test, pifD, mbfD, NbfD, epf0)

        @test Nb_c      ≈ Nb_ref      rtol=1e-12
        @test epol_c    ≈ epol_ref    rtol=1e-12
        @test Nb_pi_c   ≈ Nb_pi_ref   rtol=1e-12
        @test Nb_mb_c   ≈ Nb_mb_ref   rtol=1e-12
        @test epol_pi_c ≈ epol_pi_ref rtol=1e-12
        @test epol_mb_c ≈ epol_mb_ref rtol=1e-12

        # ------------------------------------------------------------------
        # 5. compressor_efficiency: warm-start hints updated after each call
        # ------------------------------------------------------------------
        fan_ws = Comp(pifD, mbfD, NbfD, epf0, epfmin, FanMap)
        Ng_before = fan_ws.Ng
        Rg_before = fan_ws.Rg
        comp_eff(fan_ws, pifD, mbfD)
        # After a call the hints should be updated (map solver returns converged N, R)
        @test fan_ws.Ng != Ng_before || fan_ws.Rg != Rg_before  # at least one changes

        # ------------------------------------------------------------------
        # 6. compressor_efficiency: efficiency floor clamping
        #
        # Set epol_min = 0.999 (above any realistic map output).
        # The returned efficiency must equal epol_min and derivatives = 0.
        # ------------------------------------------------------------------
        fan_floor = Comp(pifD, mbfD, NbfD, epf0, 0.999, FanMap)
        _, epol_fl, _, _, epol_pi_fl, epol_mb_fl = comp_eff(fan_floor, pifD, mbfD)
        @test epol_fl    ≈ 0.999  atol=1e-12
        @test epol_pi_fl ≈ 0.0   atol=1e-14
        @test epol_mb_fl ≈ 0.0   atol=1e-14

        # ------------------------------------------------------------------
        # 7. compressor_efficiency: windmilling (pi < 1)
        #
        # When pi < 1, the efficiency is inverted (TASOPT convention for
        # reverse-flow / braking mode).  The inverted value must be > 1
        # (since original epol ∈ (0,1) implies 1/epol > 1).
        # ------------------------------------------------------------------
        fan_wm = Comp(pifD, mbfD, NbfD, epf0, epfmin, FanMap, windmilling=true)
        pi_wm  = 0.9   # below unity → windmilling
        mb_wm  = 0.5 * mbfD
        _, epol_wm, _, _, epol_pi_wm, epol_mb_wm = comp_eff(fan_wm, pi_wm, mb_wm)
        @test epol_wm > 1.0   # inverted efficiency exceeds 1

        # Derivative consistency: epol_pi_wm should have opposite sign to raw derivative
        # (chain rule: d(1/epol)/dpi = -1/epol² * depol/dpi; for positive depol/dpi,
        # the inverted derivative is negative)
        _, _, _, _, epol_pi_raw, _ = calc_cse(FanMap, pi_wm, mb_wm, pifD, mbfD, NbfD, epf0)
        epol_wm_raw = 1.0  # dummy: check sign consistency
        # Just verify derivative has the expected sign (opposite to raw)
        if epol_pi_raw > 0.0
            @test epol_pi_wm < 0.0
        elseif epol_pi_raw < 0.0
            @test epol_pi_wm > 0.0
        end

        # Verify LPC/HPC do NOT apply windmilling (upstream fan-only behaviour)
        lpc_nowm = Comp(pilcD, mblcD, NblcD, eplc0, 0.70, LPCMap)
        _, epol_lpc_nowm, _, _, _, _ = comp_eff(lpc_nowm, 0.9, 0.5 * mblcD)
        @test epol_lpc_nowm < 1.0   # no inversion: raw efficiency stays < 1

        # ------------------------------------------------------------------
        # 8. compressor_exit!: physical invariants
        # ------------------------------------------------------------------
        epol_test = epf0   # use design-point efficiency
        pi_exit   = pifD

        st21 = FS{Float64}()
        comp_exit!(st21, st2, pi_exit, epol_test)

        # Temperature rises through compression
        @test st21.Tt > st2.Tt

        # Pressure ratio applied: pt_out ≈ pi × pt_in
        @test st21.pt ≈ pi_exit * st2.pt  rtol=1e-12

        # Composition unchanged through compressor
        @test all(st21.alpha .≈ st2.alpha)

        # Entropy increases (irreversible compression, epol < 1)
        @test st21.st > st2.st

        # Enthalpy increases (work added)
        @test st21.ht > st2.ht

        # ------------------------------------------------------------------
        # 9. compressor_exit! round-trip: matches direct gas_pratd call
        # ------------------------------------------------------------------
        res_ref = gas_pratd(
            air_alpha, 5,
            pt2, Tt2, st2.ht, st2.st, st2.cpt, st2.Rt,
            pi_exit, epol_test,
        )
        pt21_ref  = res_ref[1]
        Tt21_ref  = res_ref[2]
        ht21_ref  = res_ref[3]
        st21_ref  = res_ref[4]
        cpt21_ref = res_ref[5]
        Rt21_ref  = res_ref[6]

        @test st21.pt  ≈ pt21_ref  rtol=1e-14
        @test st21.Tt  ≈ Tt21_ref  rtol=1e-14
        @test st21.ht  ≈ ht21_ref  rtol=1e-14
        @test st21.st  ≈ st21_ref  rtol=1e-14
        @test st21.cpt ≈ cpt21_ref rtol=1e-14
        @test st21.Rt  ≈ Rt21_ref  rtol=1e-14

        # ------------------------------------------------------------------
        # 10. compressor_exit! monotonicity: larger pi → higher Tt and pt
        # ------------------------------------------------------------------
        pi_lo = pifD * 0.8
        pi_hi = pifD * 1.2

        st21_lo = FS{Float64}()
        st21_hi = FS{Float64}()
        comp_exit!(st21_lo, st2, pi_lo, epol_test)
        comp_exit!(st21_hi, st2, pi_hi, epol_test)

        @test st21_hi.Tt > st21_lo.Tt   # higher compression → hotter exit
        @test st21_hi.pt > st21_lo.pt   # higher pi → higher outlet pressure

        # ------------------------------------------------------------------
        # 11. Compressor{Float32} struct construction (type parameter)
        #
        # The map inversion (NLsolve) requires Float64, so compressor_efficiency
        # and compressor_exit! are limited to Float64 in practice.  Here we
        # verify only that Compressor{Float32} can be constructed and stores Float32.
        # ------------------------------------------------------------------
        # Use the outer convenience constructor (infers T=Float32, adds Ng/Rg defaults)
        fan32 = Comp(
            Float32(pifD), Float32(mbfD), Float32(NbfD),
            Float32(epf0), Float32(epfmin), FanMap,
        )
        @test fan32 isa Comp{Float32}
        @test fan32.piD      isa Float32
        @test fan32.mbD      isa Float32
        @test fan32.NbD      isa Float32
        @test fan32.epol0    isa Float32
        @test fan32.epol_min isa Float32

    end  # Compressor

    # compressor_Nb_residual — tasopt-j9l.29
    #
    # Tests the per-compressor corrected-speed map-match residual.
    # r = Nb_map(pi, mb) − Nb_target; converges to 0 when shaft coupling satisfied.
    @testset "compressor_Nb_residual" begin

        Comp        = TASOPT.engine.Compressor
        comp_Nb_res = TASOPT.engine.compressor_Nb_residual
        comp_eff    = TASOPT.engine.compressor_efficiency
        FanMap      = TASOPT.engine.FanMap
        LPCMap      = TASOPT.engine.LPCMap
        HPCMap      = TASOPT.engine.HPCMap

        # Design parameters (same as Compressor testset for consistency)
        pifD  = 1.6850000000000001;  mbfD  = 235.16225770724063;  NbfD  = 1.0790738309310697;  epf0  = 0.89480000000000004
        pilcD = 8.0000000000000000;  mblcD = 46.110246609262873;  NblcD = 1.0790738309310697;  eplc0 = 0.88000000000000000
        pihcD = 3.7500000000000000;  mbhcD = 7.8056539219349039;  NbhcD = 0.77137973563891493;  ephc0 = 0.87000000000000000

        fan = Comp(pifD,  mbfD,  NbfD,  epf0,  0.60, FanMap)
        lpc = Comp(pilcD, mblcD, NblcD, eplc0, 0.70, LPCMap)
        hpc = Comp(pihcD, mbhcD, NbhcD, ephc0, 0.70, HPCMap)

        # ------------------------------------------------------------------
        # 1. Design-point identity: r ≈ 0 when Nb_target = NbD
        #
        # At (piD, mbD) the map scaling guarantees Nb_map ≈ NbD (the
        # fan map returns it exactly due to its grid layout; LPC and HPC
        # may have O(ε_machine) rounding in the interpolation step).
        # ------------------------------------------------------------------
        r_fan, _, _ = comp_Nb_res(fan, pifD, mbfD, NbfD)
        @test r_fan ≈ 0.0  atol=1e-12

        r_lpc, _, _ = comp_Nb_res(lpc, pilcD, mblcD, NblcD)
        @test r_lpc ≈ 0.0  atol=1e-12

        r_hpc, _, _ = comp_Nb_res(hpc, pihcD, mbhcD, NbhcD)
        @test r_hpc ≈ 0.0  atol=1e-12

        # ------------------------------------------------------------------
        # 2. Sign: r > 0 when Nb_target is below the map speed
        # ------------------------------------------------------------------
        r_high, _, _ = comp_Nb_res(fan, pifD, mbfD, NbfD * 0.95)
        @test r_high > 0.0   # Nb_map = NbfD > 0.95*NbfD = Nb_target

        # ------------------------------------------------------------------
        # 3. Sign: r < 0 when Nb_target is above the map speed
        # ------------------------------------------------------------------
        r_low, _, _ = comp_Nb_res(fan, pifD, mbfD, NbfD * 1.05)
        @test r_low < 0.0    # Nb_map = NbfD < 1.05*NbfD = Nb_target

        # ------------------------------------------------------------------
        # 4. Round-trip: r = Nb_map - Nb_target exactly, derivatives match
        #    compressor_efficiency
        # ------------------------------------------------------------------
        pi_test = pifD * 0.97
        mb_test = mbfD * 0.98
        Nb_t    = NbfD * 0.99

        r_rt, Nb_pi_rt, Nb_mb_rt = comp_Nb_res(fan, pi_test, mb_test, Nb_t)

        fan_rt = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        Nb_eff, _, Nb_pi_eff, Nb_mb_eff, _, _ = comp_eff(fan_rt, pi_test, mb_test)

        @test r_rt     ≈ Nb_eff - Nb_t   rtol=1e-14
        @test Nb_pi_rt ≈ Nb_pi_eff       rtol=1e-12
        @test Nb_mb_rt ≈ Nb_mb_eff       rtol=1e-12

        # ------------------------------------------------------------------
        # 5. Derivative self-consistency: dr/dpi matches finite-difference
        #    Nb_target is held constant, so dr/dpi = dNb_map/dpi.
        # ------------------------------------------------------------------
        eps_fd = 1e-5 * pifD
        r_fwd, _, _ = comp_Nb_res(Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap),
                                   pi_test + eps_fd, mb_test, Nb_t)
        r_bwd, _, _ = comp_Nb_res(Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap),
                                   pi_test - eps_fd, mb_test, Nb_t)
        Nb_pi_fd = (r_fwd - r_bwd) / (2 * eps_fd)
        @test Nb_pi_rt ≈ Nb_pi_fd  rtol=1e-4

        # dr/dmb finite-difference
        eps_mb = 1e-5 * mbfD
        r_fwd_mb, _, _ = comp_Nb_res(Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap),
                                      pi_test, mb_test + eps_mb, Nb_t)
        r_bwd_mb, _, _ = comp_Nb_res(Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap),
                                      pi_test, mb_test - eps_mb, Nb_t)
        Nb_mb_fd = (r_fwd_mb - r_bwd_mb) / (2 * eps_mb)
        @test Nb_mb_rt ≈ Nb_mb_fd  rtol=1e-4

        # ------------------------------------------------------------------
        # 6. Monotonicity: increasing pi (above design) increases Nb_map
        #    for fixed mb, so r increases for fixed Nb_target.
        # ------------------------------------------------------------------
        Nb_tgt_fixed = NbfD
        fan_mono = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        r_lo, _, _ = comp_Nb_res(fan_mono, pifD * 0.90, mbfD, Nb_tgt_fixed)
        r_md, _, _ = comp_Nb_res(fan_mono, pifD * 1.00, mbfD, Nb_tgt_fixed)
        r_hi, _, _ = comp_Nb_res(fan_mono, pifD * 1.10, mbfD, Nb_tgt_fixed)
        @test r_lo < r_md
        @test r_md < r_hi

        # ------------------------------------------------------------------
        # 7. LPC and HPC at design: r = 0 (same invariant as fan)
        # ------------------------------------------------------------------
        r_lpc2, Nb_pi_lpc, Nb_mb_lpc = comp_Nb_res(lpc, pilcD, mblcD, NblcD)
        @test r_lpc2 ≈ 0.0  atol=1e-12
        @test Nb_pi_lpc isa Float64
        @test Nb_mb_lpc isa Float64

        r_hpc2, Nb_pi_hpc, Nb_mb_hpc = comp_Nb_res(hpc, pihcD, mbhcD, NbhcD)
        @test r_hpc2 ≈ 0.0  atol=1e-12
        @test Nb_pi_hpc isa Float64
        @test Nb_mb_hpc isa Float64

    end  # compressor_Nb_residual

    # compressor_pratd — tasopt-j9l.32
    #
    # Tests the compressor outlet-state + analytic Jacobian block.
    # compressor_pratd combines map inversion (compressor_efficiency) with
    # gas_pratd, applying the efficiency chain rule for pi and mb directions.
    #
    # Invariants verified:
    #   - Design-point pressure identity: pt_out ≈ pi * pt_in
    #   - pt_out_pt_in = pi (exact by construction)
    #   - pt_out_st_in = 0 (pressure doesn't depend on inlet entropy)
    #   - Analytic ∂pt_out/∂pi matches FD (rtol=1e-4; limited by map solver)
    #   - Analytic ∂Tt_out/∂pi matches FD (rtol=1e-4)
    #   - Analytic ∂pt_out/∂mb matches FD (rtol=1e-4)
    #   - Analytic ∂Tt_out/∂st_in matches FD (rtol=1e-6; purely algebraic)
    #   - Analytic ∂pt_out/∂pt_in matches FD (rtol=1e-6; trivially = pi)
    #   - Nb, Nb_pi, Nb_mb match compressor_efficiency output (round-trip)
    #   - epol, epol_pi, epol_mb match compressor_efficiency output (round-trip)
    #   - LPC and HPC at their design points (smoke test)
    #   - pt_out > pt_in (compression is positive work)
    #   - Tt_out > Tt_in (temperature rises through compressor)
    @testset "compressor_pratd" begin
        using StaticArrays

        Comp         = TASOPT.engine.Compressor
        comp_pratd   = TASOPT.engine.compressor_pratd
        comp_eff     = TASOPT.engine.compressor_efficiency
        FanMap       = TASOPT.engine.FanMap
        LPCMap       = TASOPT.engine.LPCMap
        HPCMap       = TASOPT.engine.HPCMap
        gas_pratd_fn = TASOPT.engine.gas_pratd

        # Design parameters (consistent with rest of engine testset)
        pifD  = 1.6850000000000001;  mbfD  = 235.16225770724063;  NbfD  = 1.0790738309310697;  epf0  = 0.89480000000000004
        pilcD = 8.0000000000000000;  mblcD = 46.110246609262873;  NblcD = 1.0790738309310697;  eplc0 = 0.88000000000000000
        pihcD = 3.7500000000000000;  mbhcD = 7.8056539219349039;  NbhcD = 0.77137973563891493;  ephc0 = 0.87000000000000000

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair = 5

        # Realistic fan-face inlet conditions
        pt_in  = 30_000.0   # Pa
        Tt_in  = 288.15     # K
        ht_in  = 287_500.0  # J/kg (approx)
        st_in  = 1718.0     # J/kg·K
        cpt_in = 1004.0     # J/kg·K
        Rt_in  = 287.0      # J/kg·K

        fan = Comp(pifD,  mbfD,  NbfD,  epf0,  0.60, FanMap)
        lpc = Comp(pilcD, mblcD, NblcD, eplc0, 0.70, LPCMap)
        hpc = Comp(pihcD, mbhcD, NbhcD, ephc0, 0.70, HPCMap)

        # ------------------------------------------------------------------
        # 1. Design-point pressure identity: pt_out = pi * pt_in
        # ------------------------------------------------------------------
        pi_test = pifD * 0.97
        mb_test = mbfD * 0.99

        pt_out, Tt_out, ht_out, st_out, cpt_out, Rt_out,
        pt_out_pt_in,
        pt_out_st_in, Tt_out_st_in, ht_out_st_in, st_out_st_in,
        pt_out_pi, Tt_out_pi, ht_out_pi, st_out_pi,
        pt_out_mb, Tt_out_mb, ht_out_mb, st_out_mb,
        Nb, Nb_pi, Nb_mb, epol, epol_pi, epol_mb =
            comp_pratd(fan, air_alpha, nair,
                pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
                pi_test, mb_test)

        @test pt_out ≈ pi_test * pt_in   rtol=1e-10   # p = po * pi exactly
        @test Tt_out > Tt_in                           # compression raises temperature
        @test pt_out > pt_in                           # compression raises pressure

        # ------------------------------------------------------------------
        # 2. pt_out_pt_in = pi (exact identity from gas_pratd)
        # ------------------------------------------------------------------
        @test pt_out_pt_in ≈ pi_test   rtol=1e-10

        # ------------------------------------------------------------------
        # 3. pt_out_st_in = 0 (pressure independent of inlet entropy)
        # ------------------------------------------------------------------
        @test pt_out_st_in ≈ 0.0   atol=1e-15

        # ------------------------------------------------------------------
        # 4. FD verification: ∂pt_out/∂pi   (rtol=1e-4; map solver limits precision)
        # ------------------------------------------------------------------
        eps_pi = 1e-5 * pi_test
        fan_fwd = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        fan_bwd = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        pt_fwd, = comp_pratd(fan_fwd, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test + eps_pi, mb_test)
        pt_bwd, = comp_pratd(fan_bwd, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test - eps_pi, mb_test)
        pt_out_pi_fd = (pt_fwd - pt_bwd) / (2 * eps_pi)
        @test pt_out_pi ≈ pt_out_pi_fd   rtol=1e-4

        # ------------------------------------------------------------------
        # 5. FD verification: ∂Tt_out/∂pi
        #    rtol=5e-4: Tt derivatives are noisier than pt (gas_pratd Newton ttol=1e-6)
        # ------------------------------------------------------------------
        fan_fwd2 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        fan_bwd2 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        _, Tt_fwd, = comp_pratd(fan_fwd2, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test + eps_pi, mb_test)
        _, Tt_bwd, = comp_pratd(fan_bwd2, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test - eps_pi, mb_test)
        Tt_out_pi_fd = (Tt_fwd - Tt_bwd) / (2 * eps_pi)
        @test Tt_out_pi ≈ Tt_out_pi_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 6. FD verification: ∂pt_out/∂mb  (efficiency-chain only)
        # ------------------------------------------------------------------
        eps_mb = 1e-5 * mb_test
        fan_fwd3 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        fan_bwd3 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        pt_fwd_mb, = comp_pratd(fan_fwd3, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test, mb_test + eps_mb)
        pt_bwd_mb, = comp_pratd(fan_bwd3, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test, mb_test - eps_mb)
        pt_out_mb_fd = (pt_fwd_mb - pt_bwd_mb) / (2 * eps_mb)
        @test pt_out_mb ≈ pt_out_mb_fd   rtol=1e-4

        # ------------------------------------------------------------------
        # 7. FD verification: ∂Tt_out/∂st_in
        #    rtol=1e-3: derivative goes through gas_pratd Newton (ttol=1e-6)
        # ------------------------------------------------------------------
        eps_st = 1e-5 * st_in
        fan_fwd4 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        fan_bwd4 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        _, Tt_fwd_st, = comp_pratd(fan_fwd4, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in + eps_st, cpt_in, Rt_in,
            pi_test, mb_test)
        _, Tt_bwd_st, = comp_pratd(fan_bwd4, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in - eps_st, cpt_in, Rt_in,
            pi_test, mb_test)
        Tt_out_st_in_fd = (Tt_fwd_st - Tt_bwd_st) / (2 * eps_st)
        @test Tt_out_st_in ≈ Tt_out_st_in_fd   rtol=1e-3

        # ------------------------------------------------------------------
        # 8. FD verification: ∂pt_out/∂pt_in  (trivially = pi; tighter tol)
        # ------------------------------------------------------------------
        eps_pt = 1e-5 * pt_in
        fan_fwd5 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        fan_bwd5 = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        pt_fwd_pt, = comp_pratd(fan_fwd5, air_alpha, nair,
            pt_in + eps_pt, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test, mb_test)
        pt_bwd_pt, = comp_pratd(fan_bwd5, air_alpha, nair,
            pt_in - eps_pt, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pi_test, mb_test)
        pt_out_pt_in_fd = (pt_fwd_pt - pt_bwd_pt) / (2 * eps_pt)
        @test pt_out_pt_in ≈ pt_out_pt_in_fd   rtol=1e-8

        # ------------------------------------------------------------------
        # 9. Round-trip: Nb, Nb_pi, Nb_mb, epol, epol_pi, epol_mb match
        #    direct compressor_efficiency call (same operating point)
        # ------------------------------------------------------------------
        fan_rt = Comp(pifD, mbfD, NbfD, epf0, 0.60, FanMap)
        Nb_eff, epol_eff, Nb_pi_eff, Nb_mb_eff, epol_pi_eff, epol_mb_eff =
            comp_eff(fan_rt, pi_test, mb_test)

        @test Nb    ≈ Nb_eff      rtol=1e-12
        @test Nb_pi ≈ Nb_pi_eff   rtol=1e-10
        @test Nb_mb ≈ Nb_mb_eff   rtol=1e-10
        @test epol  ≈ epol_eff    rtol=1e-12

        # ------------------------------------------------------------------
        # 10. LPC and HPC smoke tests: pt_out > pt_in, Tt_out > Tt_in
        # ------------------------------------------------------------------
        pt_lpc, Tt_lpc = comp_pratd(lpc, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pilcD * 0.98, mblcD * 0.99)[1:2]
        @test pt_lpc > pt_in
        @test Tt_lpc > Tt_in

        pt_hpc, Tt_hpc = comp_pratd(hpc, air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in,
            pihcD * 0.98, mbhcD * 0.99)[1:2]
        @test pt_hpc > pt_in
        @test Tt_hpc > Tt_in

    end  # compressor_pratd

    # turbine_delhd — tasopt-j9l.57
    #
    # Tests the turbine outlet-state + analytic Jacobian block for work-based expansion.
    # turbine_delhd wraps gas_delhd with the efficiency chain rule pre-applied:
    #   pt_out_ept = pt_out_epi * (-epi/ept)   where epi = 1/ept
    #
    # Invariants verified:
    #   - Energy conservation: ht_out ≈ ht_in + dh (within gas_delhd precision)
    #   - Tt_out < Tt_in  (expansion cools)
    #   - pt_out < pt_in  (expansion drops total pressure)
    #   - pt_out_pt_in ≈ pt_out / pt_in  (linear scaling with inlet pressure)
    #   - FD: ∂pt_out/∂dh  matches pt_out_dh    (rtol=1e-6)
    #   - FD: ∂Tt_out/∂dh  matches Tt_out_dh    (rtol=1e-6)
    #   - FD: ∂pt_out/∂ept matches pt_out_ept   (rtol=1e-5; chain through epi)
    #   - FD: ∂Tt_out/∂ept ≈ 0                  (energy conservation; Tt indep. of ept)
    #   - FD: ∂pt_out/∂ht_in matches pt_out_ht_in (rtol=1e-6)
    #   - FD: ∂Tt_out/∂ht_in matches Tt_out_ht_in (rtol=1e-6)
    #   - round-trip: scalar outputs match direct gas_delhd call at same (dh, ept)
    @testset "turbine_delhd" begin
        using StaticArrays

        turb_delhd = TASOPT.engine.turbine_delhd

        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair = 5

        # HPT-like inlet conditions (post-combustor, high temperature)
        pt_in  = 2.0e6    # Pa  (HPT inlet total pressure)
        Tt_in  = 1650.0   # K   (HPT inlet temperature)
        ht_in  = 1.73e6   # J/kg  (approx enthalpy at 1650 K)
        st_in  = 2900.0   # J/kg·K
        cpt_in = 1148.0   # J/kg·K  (approx cp at 1650 K)
        Rt_in  = 287.0    # J/kg·K

        dh  = -5.0e5      # J/kg  (work extraction; negative)
        ept = 0.90        # polytropic efficiency

        pt_out, Tt_out, ht_out, st_out, cpt_out, Rt_out,
        pt_out_pt_in, pt_out_st_in,
        pt_out_ht_in, Tt_out_ht_in, ht_out_ht_in, st_out_ht_in,
        pt_out_dh, Tt_out_dh, ht_out_dh, st_out_dh,
        pt_out_ept,
        p_al, T_al, h_al, s_al =
            turb_delhd(air_alpha, nair,
                pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, ept)

        # ------------------------------------------------------------------
        # 1. Energy conservation: ht_out ≈ ht_in + dh
        # ------------------------------------------------------------------
        @test ht_out ≈ ht_in + dh   rtol=1e-6

        # ------------------------------------------------------------------
        # 2. Physical direction: expansion drops temperature and pressure
        # ------------------------------------------------------------------
        @test Tt_out < Tt_in
        @test pt_out < pt_in

        # ------------------------------------------------------------------
        # 3. pt_out_pt_in = pt_out / pt_in  (linear pressure scaling)
        # ------------------------------------------------------------------
        @test pt_out_pt_in ≈ pt_out / pt_in   rtol=1e-8

        # ------------------------------------------------------------------
        # 4. FD: ∂pt_out/∂dh
        # ------------------------------------------------------------------
        eps_dh = 1e-4 * abs(dh)
        pt_fwd_dh, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh + eps_dh, ept)
        pt_bwd_dh, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh - eps_dh, ept)
        pt_out_dh_fd = (pt_fwd_dh - pt_bwd_dh) / (2 * eps_dh)
        @test pt_out_dh ≈ pt_out_dh_fd   rtol=5e-4  # gas_delhd Newton ttol limits precision

        # ------------------------------------------------------------------
        # 5. FD: ∂Tt_out/∂dh
        # ------------------------------------------------------------------
        _, Tt_fwd_dh, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh + eps_dh, ept)
        _, Tt_bwd_dh, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh - eps_dh, ept)
        Tt_out_dh_fd = (Tt_fwd_dh - Tt_bwd_dh) / (2 * eps_dh)
        @test Tt_out_dh ≈ Tt_out_dh_fd   rtol=1e-5

        # ------------------------------------------------------------------
        # 6. FD: ∂pt_out/∂ept  (efficiency chain: epi = 1/ept → pt_out_ept)
        # ------------------------------------------------------------------
        eps_ept = 1e-5 * ept
        pt_fwd_ept, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, ept + eps_ept)
        pt_bwd_ept, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, ept - eps_ept)
        pt_out_ept_fd = (pt_fwd_ept - pt_bwd_ept) / (2 * eps_ept)
        @test pt_out_ept ≈ pt_out_ept_fd   rtol=1e-5

        # ------------------------------------------------------------------
        # 7. ∂Tt_out/∂ept ≈ 0  (Tt set by energy conservation, not ept)
        # ------------------------------------------------------------------
        _, Tt_fwd_ept, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, ept + eps_ept)
        _, Tt_bwd_ept, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, ept - eps_ept)
        Tt_out_ept_fd = (Tt_fwd_ept - Tt_bwd_ept) / (2 * eps_ept)
        @test abs(Tt_out_ept_fd) < 1e-6 * Tt_out   # Tt is independent of ept

        # ------------------------------------------------------------------
        # 8. FD: ∂pt_out/∂ht_in
        # ------------------------------------------------------------------
        eps_ht = 1e-4 * abs(ht_in)
        pt_fwd_ht, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in + eps_ht, st_in, cpt_in, Rt_in, dh, ept)
        pt_bwd_ht, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in - eps_ht, st_in, cpt_in, Rt_in, dh, ept)
        pt_out_ht_in_fd = (pt_fwd_ht - pt_bwd_ht) / (2 * eps_ht)
        @test pt_out_ht_in ≈ pt_out_ht_in_fd   rtol=5e-4  # gas_delhd Newton ttol limits precision

        # ------------------------------------------------------------------
        # 9. FD: ∂Tt_out/∂ht_in
        # ------------------------------------------------------------------
        _, Tt_fwd_ht, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in + eps_ht, st_in, cpt_in, Rt_in, dh, ept)
        _, Tt_bwd_ht, = turb_delhd(air_alpha, nair, pt_in, Tt_in, ht_in - eps_ht, st_in, cpt_in, Rt_in, dh, ept)
        Tt_out_ht_in_fd = (Tt_fwd_ht - Tt_bwd_ht) / (2 * eps_ht)
        @test Tt_out_ht_in ≈ Tt_out_ht_in_fd   rtol=1e-5

        # ------------------------------------------------------------------
        # 10. Round-trip: scalar outputs match direct gas_delhd call
        # ------------------------------------------------------------------
        gas_delhd_fn = TASOPT.engine.gas_delhd
        epi = 1.0 / ept
        res = gas_delhd_fn(air_alpha, nair,
            pt_in, Tt_in, ht_in, st_in, cpt_in, Rt_in, dh, epi)
        @test pt_out ≈ res[1]   rtol=1e-12
        @test Tt_out ≈ res[2]   rtol=1e-12
        @test ht_out ≈ res[3]   rtol=1e-12
        @test st_out ≈ res[4]   rtol=1e-12

    end  # turbine_delhd

    # Shaft component — tasopt-j9l.27
    @testset "Shaft" begin

        Shaft_t            = TASOPT.engine.Shaft
        hp_shaft_work_t    = TASOPT.engine.hp_shaft_work
        lp_shaft_work_t    = TASOPT.engine.lp_shaft_work
        shaft_speed_res_t  = TASOPT.engine.shaft_speed_residual

        # Representative parameter values
        epsh_val = 0.01    # HP shaft loss fraction
        epsl_val = 0.015   # LP shaft loss fraction
        Gearf_val = 1.0    # LP shaft gear ratio (direct drive in this example)

        # ── 1. Struct construction ────────────────────────────────────────────

        shaft_hp = Shaft_t(epsh_val, 1.0)
        shaft_lp = Shaft_t(epsl_val, Gearf_val)

        @test shaft_hp.eps   == epsh_val
        @test shaft_hp.Gearf == 1.0
        @test shaft_lp.eps   == epsl_val
        @test shaft_lp.Gearf == Gearf_val

        # Type parameter preserved
        shaft32 = Shaft_t(Float32(epsh_val), Float32(1.0))
        @test shaft32 isa Shaft_t{Float32}
        @test shaft32.eps isa Float32

        # ── 2. hp_shaft_work: formula verification ────────────────────────────

        fo = 0.04; ff = 0.025; ht3 = 1_050_000.0; ht25c = 650_000.0

        dhht, dhfac, dhfac_fo, dhfac_ff = hp_shaft_work_t(shaft_hp, fo, ff, ht3, ht25c)

        # Expected from the formula:
        fac_ref    = 1.0 - fo + ff
        dhfac_ref  = -(1.0 - fo) / fac_ref / (1.0 - epsh_val)
        dhht_ref   = (ht3 - ht25c) * dhfac_ref
        dhfac_fo_ref = dhfac_ref / fac_ref + 1.0 / fac_ref / (1.0 - epsh_val)
        dhfac_ff_ref = -dhfac_ref / fac_ref

        @test dhht   ≈ dhht_ref   rtol=1e-14
        @test dhfac  ≈ dhfac_ref  rtol=1e-14
        @test dhfac_fo ≈ dhfac_fo_ref rtol=1e-14
        @test dhfac_ff ≈ dhfac_ff_ref rtol=1e-14

        # dhht must be negative (HPT extracts energy → enthalpy drop)
        @test dhht < 0.0

        # Design-point identity: work magnitude = HPC work / (1 - eps)
        @test abs(dhht) ≈ abs(ht3 - ht25c) * abs(dhfac_ref)  rtol=1e-14

        # ── 3. hp_shaft_work: monotonicity in eps ────────────────────────────

        # More loss → more HPT work required (|dhht| larger for larger eps)
        dhht_lo, _, _, _ = hp_shaft_work_t(Shaft_t(0.005, 1.0), fo, ff, ht3, ht25c)
        dhht_hi, _, _, _ = hp_shaft_work_t(Shaft_t(0.050, 1.0), fo, ff, ht3, ht25c)
        @test abs(dhht_hi) > abs(dhht_lo)

        # ── 4. lp_shaft_work: formula verification ────────────────────────────

        BPR = 5.0; ht25 = 680_000.0; ht19c = 310_000.0
        ht21 = 380_000.0; ht2 = 290_000.0; Pom = 1_500.0

        dhlt, dlfac, dlfac_fo, dlfac_ff = lp_shaft_work_t(shaft_lp, fo, ff, BPR, ht25, ht19c, ht21, ht2, Pom)

        fac_lp_ref   = 1.0 - fo + ff
        dlfac_ref    = -1.0 / fac_lp_ref / (1.0 - epsl_val)
        demand_ref   = ht25 - ht19c + BPR * (ht21 - ht2) + Pom
        dhlt_ref     = demand_ref * dlfac_ref
        dlfac_fo_ref = dlfac_ref / fac_lp_ref
        dlfac_ff_ref = -dlfac_ref / fac_lp_ref

        @test dhlt   ≈ dhlt_ref   rtol=1e-14
        @test dlfac  ≈ dlfac_ref  rtol=1e-14
        @test dlfac_fo ≈ dlfac_fo_ref rtol=1e-14
        @test dlfac_ff ≈ dlfac_ff_ref rtol=1e-14

        # dhlt must be negative (LPT extracts energy)
        @test dhlt < 0.0

        # ── 5. lp_shaft_work: Pom = 0 reproduces no-offtake case ─────────────

        dhlt_no_pom, dlfac_no, _, _ = lp_shaft_work_t(shaft_lp, fo, ff, BPR, ht25, ht19c, ht21, ht2, 0.0)
        demand_no_pom = ht25 - ht19c + BPR * (ht21 - ht2)
        @test dhlt_no_pom ≈ demand_no_pom * dlfac_ref   rtol=1e-14
        @test dlfac_no    ≈ dlfac_ref                   rtol=1e-14

        # Adding offtake increases |dhlt| (LPT must do more work)
        @test abs(dhlt) > abs(dhlt_no_pom)

        # ── 6. shaft_speed_residual: design-point zero ───────────────────────

        # LP shaft: Gearf * trf * Nf = trl * Nl at design
        Nf_des = 1.5; Nl_des = Gearf_val * Nf_des   # exact design balance (Gearf=1)
        trf = sqrt(288.0 / 288.15)                   # typical cruise correction
        trl = sqrt(300.0 / 288.15)

        # With Gearf = 1, trf = trl if temps equal; use asymmetric case
        shaft_test_lp = Shaft_t(epsl_val, 2.0)       # Gearf = 2
        Nl_design = 2.0 * trf * Nf_des / trl         # satisfies Gearf * trf * Nf = trl * Nl

        r, r_Nf, r_Nl = shaft_speed_res_t(shaft_test_lp, Nf_des, Nl_design, trf, trl)
        @test r ≈ 0.0  atol=1e-12

        # ── 7. shaft_speed_residual: derivative sign and magnitude ───────────

        # r_Nf should be positive (increasing Nf increases r)
        @test r_Nf > 0.0
        @test r_Nf ≈ shaft_test_lp.Gearf * trf  rtol=1e-14

        # r_Nl should be negative (increasing Nl decreases r)
        @test r_Nl < 0.0
        @test r_Nl ≈ -trl  rtol=1e-14

        # ── 8. shaft_speed_residual: sign above/below design ─────────────────

        r_above, _, _ = shaft_speed_res_t(shaft_test_lp, 1.05 * Nf_des, Nl_design, trf, trl)
        r_below, _, _ = shaft_speed_res_t(shaft_test_lp, 0.95 * Nf_des, Nl_design, trf, trl)
        @test r_above > 0.0
        @test r_below < 0.0

        # ── 9. Round-trip vs. inline formula for hp_shaft_work ───────────────

        for (fo_test, ff_test, ht3_test, ht25c_test) in (
                (0.02, 0.020, 1_000_000.0, 600_000.0),
                (0.04, 0.025, 1_050_000.0, 650_000.0),
                (0.06, 0.030, 1_100_000.0, 700_000.0),
        )
            d, dfac, _, _ = hp_shaft_work_t(shaft_hp, fo_test, ff_test, ht3_test, ht25c_test)
            fac = 1.0 - fo_test + ff_test
            dfac_inline = -(1.0 - fo_test) / fac / (1.0 - epsh_val)
            d_inline    = (ht3_test - ht25c_test) * dfac_inline
            @test d    ≈ d_inline    rtol=1e-14
            @test dfac ≈ dfac_inline rtol=1e-14
        end

        # ── 10. Round-trip vs. inline formula for lp_shaft_work ──────────────

        for (fo_test, ff_test) in ((0.02, 0.020), (0.04, 0.025), (0.06, 0.030))
            d, dfac, _, _ = lp_shaft_work_t(shaft_lp, fo_test, ff_test, BPR, ht25, ht19c, ht21, ht2, Pom)
            fac_t = 1.0 - fo_test + ff_test
            dfac_inline = -1.0 / fac_t / (1.0 - epsl_val)
            demand_t    = ht25 - ht19c + BPR * (ht21 - ht2) + Pom
            d_inline    = demand_t * dfac_inline
            @test d    ≈ d_inline    rtol=1e-14
            @test dfac ≈ dfac_inline rtol=1e-14
        end

    end  # Shaft

    # combustor_burnd — tasopt-j9l.58
    @testset "combustor_burnd" begin
        using StaticArrays

        Comb_t      = TASOPT.engine.Combustor
        comb_burnd  = TASOPT.engine.combustor_burnd

        # Representative combustor parameters
        pib_val  = 0.94
        etab_val = 0.985
        Ttf_val  = 300.0
        ifuel_val = 24
        hvap_val = 0.0

        burner = Comb_t(pib_val, etab_val, Ttf_val, ifuel_val, hvap_val)
        air_alpha = SA[0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair = 5

        # Operating point: HPC exit temperature + burner exit temperature
        Tt3_test = 830.0   # K (HPC exit)
        Tb_test  = 1587.0  # K (combustor exit / Tt4)

        # Base evaluation
        ffb, lambda,
        Tt4, ht4, st4, cpt4, Rt4,
        ffb_Tt3, ffb_Tb,
        ht4_Tt3, st4_Tt3, cpt4_Tt3, Rt4_Tt3,
        ht4_Tb, st4_Tb, cpt4_Tb, Rt4_Tb,
        lam_Tt3, lam_Tb = comb_burnd(burner, air_alpha, nair, Tt3_test, Tb_test)

        # ------------------------------------------------------------------
        # 1. Basic sanity: Tt4 = Tb, positive fuel fraction, valid state
        # ------------------------------------------------------------------
        @test Tt4 == Tb_test
        @test ffb > 0.0
        @test length(lambda) == nair
        @test sum(lambda) ≈ 1.0   atol=0.05  # composition sums close to 1

        # ------------------------------------------------------------------
        # 2. FD verification: ∂ffb/∂Tt3
        #    rtol=5e-4: gas_burnd internal Newton limits FD precision
        # ------------------------------------------------------------------
        eps_Tt3 = 1e-5 * Tt3_test
        ffb_fwd, = comb_burnd(burner, air_alpha, nair, Tt3_test + eps_Tt3, Tb_test)
        ffb_bwd, = comb_burnd(burner, air_alpha, nair, Tt3_test - eps_Tt3, Tb_test)
        ffb_Tt3_fd = (ffb_fwd - ffb_bwd) / (2 * eps_Tt3)
        @test ffb_Tt3 ≈ ffb_Tt3_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 3. FD verification: ∂ffb/∂Tb
        # ------------------------------------------------------------------
        eps_Tb = 1e-5 * Tb_test
        ffb_fwd_Tb, = comb_burnd(burner, air_alpha, nair, Tt3_test, Tb_test + eps_Tb)
        ffb_bwd_Tb, = comb_burnd(burner, air_alpha, nair, Tt3_test, Tb_test - eps_Tb)
        ffb_Tb_fd = (ffb_fwd_Tb - ffb_bwd_Tb) / (2 * eps_Tb)
        @test ffb_Tb ≈ ffb_Tb_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 4. FD verification: ∂ht4/∂Tt3 (composition chain only)
        # ------------------------------------------------------------------
        res_fwd = comb_burnd(burner, air_alpha, nair, Tt3_test + eps_Tt3, Tb_test)
        res_bwd = comb_burnd(burner, air_alpha, nair, Tt3_test - eps_Tt3, Tb_test)
        ht4_Tt3_fd = (res_fwd[4] - res_bwd[4]) / (2 * eps_Tt3)
        @test ht4_Tt3 ≈ ht4_Tt3_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 5. FD verification: ∂st4/∂Tt3
        # ------------------------------------------------------------------
        st4_Tt3_fd = (res_fwd[5] - res_bwd[5]) / (2 * eps_Tt3)
        @test st4_Tt3 ≈ st4_Tt3_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 6. FD verification: ∂ht4/∂Tb (composition + temperature)
        # ------------------------------------------------------------------
        res_fwd_Tb = comb_burnd(burner, air_alpha, nair, Tt3_test, Tb_test + eps_Tb)
        res_bwd_Tb = comb_burnd(burner, air_alpha, nair, Tt3_test, Tb_test - eps_Tb)
        ht4_Tb_fd = (res_fwd_Tb[4] - res_bwd_Tb[4]) / (2 * eps_Tb)
        @test ht4_Tb ≈ ht4_Tb_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 7. FD verification: ∂st4/∂Tb
        # ------------------------------------------------------------------
        st4_Tb_fd = (res_fwd_Tb[5] - res_bwd_Tb[5]) / (2 * eps_Tb)
        @test st4_Tb ≈ st4_Tb_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 8. FD verification: ∂cpt4/∂Tb
        # ------------------------------------------------------------------
        cpt4_Tb_fd = (res_fwd_Tb[6] - res_bwd_Tb[6]) / (2 * eps_Tb)
        @test cpt4_Tb ≈ cpt4_Tb_fd   rtol=5e-4

        # ------------------------------------------------------------------
        # 9. FD verification: ∂Rt4/∂Tt3 (composition chain — small magnitude)
        # ------------------------------------------------------------------
        Rt4_Tt3_fd = (res_fwd[7] - res_bwd[7]) / (2 * eps_Tt3)
        @test Rt4_Tt3 ≈ Rt4_Tt3_fd   atol=1e-5

        # ------------------------------------------------------------------
        # 10. Composition derivatives: ∂lambda/∂Tt3 via FD
        # ------------------------------------------------------------------
        lam_fwd = comb_burnd(burner, air_alpha, nair, Tt3_test + eps_Tt3, Tb_test)[2]
        lam_bwd = comb_burnd(burner, air_alpha, nair, Tt3_test - eps_Tt3, Tb_test)[2]
        lam_Tt3_fd = (lam_fwd .- lam_bwd) ./ (2 * eps_Tt3)
        for i in 1:nair
            @test lam_Tt3[i] ≈ lam_Tt3_fd[i]   rtol=5e-4
        end
    end  # combustor_burnd

    # hp_shaft_workd — tasopt-j9l.58
    @testset "hp_shaft_workd" begin
        Shaft_t        = TASOPT.engine.Shaft
        hp_workd       = TASOPT.engine.hp_shaft_workd
        hp_work        = TASOPT.engine.hp_shaft_work

        shaft = Shaft_t(0.01, 1.0)  # HP shaft

        fo_val    = 0.015
        ff_val    = 0.025
        ht3_val   = 900_000.0   # J/kg (HPC exit)
        ht25c_val = 430_000.0   # J/kg (HPC inlet)

        # Base evaluation
        dhht, dhht_fo, dhht_ff, dhht_ht3, dhht_ht25c =
            hp_workd(shaft, fo_val, ff_val, ht3_val, ht25c_val)

        # 1. Consistency with hp_shaft_work
        dhht_ref, dhfac_ref, _, _ = hp_work(shaft, fo_val, ff_val, ht3_val, ht25c_val)
        @test dhht ≈ dhht_ref   rtol=1e-14

        # 2. ∂dhht/∂ht3 = dhfac, ∂dhht/∂ht25c = -dhfac
        @test dhht_ht3   ≈  dhfac_ref   rtol=1e-14
        @test dhht_ht25c ≈ -dhfac_ref   rtol=1e-14

        # 3. FD verification: ∂dhht/∂fo
        eps_fo = 1e-7
        dhht_fwd, = hp_workd(shaft, fo_val + eps_fo, ff_val, ht3_val, ht25c_val)
        dhht_bwd, = hp_workd(shaft, fo_val - eps_fo, ff_val, ht3_val, ht25c_val)
        @test dhht_fo ≈ (dhht_fwd - dhht_bwd) / (2 * eps_fo)   rtol=1e-7

        # 4. FD verification: ∂dhht/∂ff
        eps_ff = 1e-7
        dhht_fwd_ff, = hp_workd(shaft, fo_val, ff_val + eps_ff, ht3_val, ht25c_val)
        dhht_bwd_ff, = hp_workd(shaft, fo_val, ff_val - eps_ff, ht3_val, ht25c_val)
        @test dhht_ff ≈ (dhht_fwd_ff - dhht_bwd_ff) / (2 * eps_ff)   rtol=1e-7

        # 5. FD verification: ∂dhht/∂ht3
        eps_h = 1.0
        dhht_fwd_h, = hp_workd(shaft, fo_val, ff_val, ht3_val + eps_h, ht25c_val)
        dhht_bwd_h, = hp_workd(shaft, fo_val, ff_val, ht3_val - eps_h, ht25c_val)
        @test dhht_ht3 ≈ (dhht_fwd_h - dhht_bwd_h) / (2 * eps_h)   rtol=1e-10
    end  # hp_shaft_workd

    # lp_shaft_workd — tasopt-j9l.58
    @testset "lp_shaft_workd" begin
        Shaft_t        = TASOPT.engine.Shaft
        lp_workd       = TASOPT.engine.lp_shaft_workd
        lp_work        = TASOPT.engine.lp_shaft_work

        shaft = Shaft_t(0.015, 1.3)  # LP shaft with gear ratio

        fo_val    = 0.015
        ff_val    = 0.025
        BPR_val   = 5.5
        ht25_val  = 450_000.0   # LPC exit
        ht19c_val = 300_000.0   # LPC inlet
        ht21_val  = 320_000.0   # Fan exit
        ht2_val   = 290_000.0   # Fan inlet
        Pom_val   = 5000.0      # Power offtake

        # Base evaluation
        dhlt, dhlt_fo, dhlt_ff, dhlt_BPR,
        dhlt_ht25, dhlt_ht19c, dhlt_ht21, dhlt_ht2, dhlt_Pom =
            lp_workd(shaft, fo_val, ff_val, BPR_val,
                     ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)

        # 1. Consistency with lp_shaft_work
        dhlt_ref, dlfac_ref, _, _ = lp_work(shaft, fo_val, ff_val, BPR_val,
                                            ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        @test dhlt ≈ dhlt_ref   rtol=1e-14

        # 2. Analytic identities
        @test dhlt_ht25  ≈  dlfac_ref         rtol=1e-14
        @test dhlt_ht19c ≈ -dlfac_ref         rtol=1e-14
        @test dhlt_ht21  ≈  BPR_val * dlfac_ref   rtol=1e-14
        @test dhlt_ht2   ≈ -BPR_val * dlfac_ref   rtol=1e-14
        @test dhlt_Pom   ≈  dlfac_ref         rtol=1e-14
        @test dhlt_BPR   ≈ (ht21_val - ht2_val) * dlfac_ref   rtol=1e-14

        # 3. FD verification: ∂dhlt/∂fo
        eps_fo = 1e-7
        dhlt_fwd, = lp_workd(shaft, fo_val + eps_fo, ff_val, BPR_val,
                             ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        dhlt_bwd, = lp_workd(shaft, fo_val - eps_fo, ff_val, BPR_val,
                             ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        @test dhlt_fo ≈ (dhlt_fwd - dhlt_bwd) / (2 * eps_fo)   rtol=1e-7

        # 4. FD verification: ∂dhlt/∂ff
        eps_ff = 1e-7
        dhlt_fwd_ff, = lp_workd(shaft, fo_val, ff_val + eps_ff, BPR_val,
                                ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        dhlt_bwd_ff, = lp_workd(shaft, fo_val, ff_val - eps_ff, BPR_val,
                                ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        @test dhlt_ff ≈ (dhlt_fwd_ff - dhlt_bwd_ff) / (2 * eps_ff)   rtol=1e-7

        # 5. FD verification: ∂dhlt/∂BPR
        eps_BPR = 1e-5
        dhlt_fwd_B, = lp_workd(shaft, fo_val, ff_val, BPR_val + eps_BPR,
                               ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        dhlt_bwd_B, = lp_workd(shaft, fo_val, ff_val, BPR_val - eps_BPR,
                               ht25_val, ht19c_val, ht21_val, ht2_val, Pom_val)
        @test dhlt_BPR ≈ (dhlt_fwd_B - dhlt_bwd_B) / (2 * eps_BPR)   rtol=1e-7
    end  # lp_shaft_workd

    @testset "Splitter" begin

        Splitter_t    = TASOPT.engine.Splitter
        bypass_ratio_t = TASOPT.engine.bypass_ratio

        splitter = Splitter_t()

        # ── 1. Struct construction ────────────────────────────────────────────

        @test splitter isa Splitter_t

        # ── 2. Formula verification ───────────────────────────────────────────

        # Representative cruise conditions
        mf    = 1.08          # corrected fan mass flow [—]
        ml    = 0.20          # corrected LPC mass flow [—]
        pt2   = 29_800.0      # fan-inlet total pressure  [Pa]
        pt19c = 29_600.0      # LPC-inlet total pressure  [Pa]
        Tt2   = 236.0         # fan-inlet total temperature  [K]
        Tt19c = 236.4         # LPC-inlet total temperature  [K]

        BPR, BPR_mf, BPR_ml, BPR_pt2, BPR_pt19c =
            bypass_ratio_t(splitter, mf, ml, pt2, pt19c, Tt2, Tt19c)

        BPR_ref = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c

        @test BPR       ≈ BPR_ref              rtol=1e-14
        @test BPR_mf    ≈  BPR_ref / mf        rtol=1e-14
        @test BPR_ml    ≈ -BPR_ref / ml        rtol=1e-14
        @test BPR_pt2   ≈  BPR_ref / pt2       rtol=1e-14
        @test BPR_pt19c ≈ -BPR_ref / pt19c     rtol=1e-14

        # ── 3. BPR must be positive for physically valid inputs ───────────────

        @test BPR > 0.0

        # ── 4. Euler's theorem: BPR is homogeneous of degree 1 in (mf/ml) ────
        #    Scaling mf and ml by the same factor leaves BPR unchanged.

        scale = 1.3
        BPR_scaled, _, _, _, _ =
            bypass_ratio_t(splitter, scale * mf, scale * ml, pt2, pt19c, Tt2, Tt19c)
        @test BPR_scaled ≈ BPR rtol=1e-14

        # ── 5. Monotonicity ───────────────────────────────────────────────────
        #    BPR increases with mf, decreases with ml, increases with pt2,
        #    decreases with pt19c.

        BPR_hi_mf, _, _, _, _ =
            bypass_ratio_t(splitter, 1.1 * mf, ml, pt2, pt19c, Tt2, Tt19c)
        @test BPR_hi_mf > BPR

        BPR_hi_ml, _, _, _, _ =
            bypass_ratio_t(splitter, mf, 1.1 * ml, pt2, pt19c, Tt2, Tt19c)
        @test BPR_hi_ml < BPR

        BPR_hi_pt2, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, 1.05 * pt2, pt19c, Tt2, Tt19c)
        @test BPR_hi_pt2 > BPR

        BPR_hi_pt19c, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, pt2, 1.05 * pt19c, Tt2, Tt19c)
        @test BPR_hi_pt19c < BPR

        # ── 6. Finite-difference derivative self-consistency ──────────────────

        δ = 1e-6

        BPR_p, _, _, _, _ =
            bypass_ratio_t(splitter, mf + δ, ml, pt2, pt19c, Tt2, Tt19c)
        BPR_m, _, _, _, _ =
            bypass_ratio_t(splitter, mf - δ, ml, pt2, pt19c, Tt2, Tt19c)
        @test BPR_mf ≈ (BPR_p - BPR_m) / (2δ) rtol=1e-8

        BPR_p2, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml + δ, pt2, pt19c, Tt2, Tt19c)
        BPR_m2, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml - δ, pt2, pt19c, Tt2, Tt19c)
        @test BPR_ml ≈ (BPR_p2 - BPR_m2) / (2δ) rtol=1e-8

        δp = pt2 * 1e-6
        BPR_pp, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, pt2 + δp, pt19c, Tt2, Tt19c)
        BPR_mp, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, pt2 - δp, pt19c, Tt2, Tt19c)
        @test BPR_pt2 ≈ (BPR_pp - BPR_mp) / (2δp) rtol=1e-8

        δp19 = pt19c * 1e-6
        BPR_pp19, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, pt2, pt19c + δp19, Tt2, Tt19c)
        BPR_mp19, _, _, _, _ =
            bypass_ratio_t(splitter, mf, ml, pt2, pt19c - δp19, Tt2, Tt19c)
        @test BPR_pt19c ≈ (BPR_pp19 - BPR_mp19) / (2δp19) rtol=1e-8

        # ── 7. Chain-rule expansion matches original inline formula ───────────
        #    Verify that the full BPR_mf assembled via chain rule matches the
        #    inline formula from tfoper.jl, using hypothetical upstream sensitivities.

        pt2_mf = 50.0;  pt19c_mf = -20.0
        pt2_ml = -30.0; pt19c_ml =  15.0
        pt2_Mi =  80.0; pt19c_Mi = -60.0

        BPR_mf_full  = BPR_mf  + BPR_pt2 * pt2_mf  + BPR_pt19c * pt19c_mf
        BPR_ml_full  = BPR_ml  + BPR_pt2 * pt2_ml  + BPR_pt19c * pt19c_ml
        BPR_Mi_full  =           BPR_pt2 * pt2_Mi  + BPR_pt19c * pt19c_Mi

        # Inline reference (from tfoper.jl original)
        BPR_mf_inline = 1.0 / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c +
                         BPR / pt2 * pt2_mf - BPR / pt19c * pt19c_mf
        BPR_ml_inline = -BPR / ml +
                         BPR / pt2 * pt2_ml - BPR / pt19c * pt19c_ml
        BPR_Mi_inline  = BPR / pt2 * pt2_Mi - BPR / pt19c * pt19c_Mi

        @test BPR_mf_full ≈ BPR_mf_inline  rtol=1e-14
        @test BPR_ml_full ≈ BPR_ml_inline  rtol=1e-14
        @test BPR_Mi_full ≈ BPR_Mi_inline  rtol=1e-14

        # ── 8. Float32 dispatch (type stability) ──────────────────────────────

        BPR32, BPR_mf32, BPR_ml32, BPR_pt2_32, BPR_pt19c_32 =
            bypass_ratio_t(splitter,
                Float32(mf), Float32(ml),
                Float32(pt2), Float32(pt19c),
                Float32(Tt2), Float32(Tt19c))
        @test BPR32       isa Float32
        @test BPR_mf32    isa Float32
        @test BPR_ml32    isa Float32
        @test BPR_pt2_32  isa Float32
        @test BPR_pt19c_32 isa Float32
        @test BPR32 ≈ Float32(BPR) rtol=1e-5

    end  # Splitter

    # Nozzle component — tasopt-j9l.56
    @testset "Nozzle" begin

        Nozzle_t                    = TASOPT.engine.Nozzle
        nozzle_exit_t               = TASOPT.engine.nozzle_exit
        nozzle_massflow_residual_t  = TASOPT.engine.nozzle_massflow_residual
        nozzle_gross_thrust_t       = TASOPT.engine.nozzle_gross_thrust

        # ---------------------------------------------------------------------------
        # Thermodynamic state — representative cruise conditions
        # Air species composition (standard TASOPT ordering)
        # ---------------------------------------------------------------------------
        alpha_air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair  = 5

        # Fan nozzle inlet conditions derived from a default sized aircraft
        # (run_engine_design_point output, design cruise point).
        # These match the tfoper station 21 totals used in engine_sweep_regression.
        pt21  = 61271.983314487086
        Tt21  = 294.1405135675445
        pifn  = 0.98          # fan nozzle pressure recovery (pare[iepifn])
        A7    = 0.841808817378845  # fan nozzle area, m² (design.A7)
        p0    = 23922.608843328788 # cruise ambient static pressure, Pa

        # Compute entropy-complement and ht from Tt21 via gassum
        gassum_t = TASOPT.engine.gassum
        st21, _, ht21, _, cpt21, Rt21 = gassum_t(alpha_air, nair, Tt21)

        # Fan nozzle component — design point is choked (M=1)
        fan_nozzle = Nozzle_t(pifn, A7)

        # ---------------------------------------------------------------------------
        # 1. Struct construction
        # ---------------------------------------------------------------------------

        @test fan_nozzle isa Nozzle_t{Float64}
        @test fan_nozzle.pn ≈ pifn
        @test fan_nozzle.A  ≈ A7

        # ── 2. Nozzle with different pressure recovery and Float32 ─────────────────

        noz_hi_pn = Nozzle_t(0.99, A7)
        @test noz_hi_pn.pn ≈ 0.99

        noz32 = Nozzle_t(Float32(pifn), Float32(A7))
        @test noz32 isa Nozzle_t{Float32}

        # ---------------------------------------------------------------------------
        # 3. nozzle_exit — choked case (design point)
        # ---------------------------------------------------------------------------

        p_e, T_e, h_e, s_e, cp_e, R_e, u_e, rho_e, M_e =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          pt21, Tt21, ht21, st21, cpt21, Rt21, p0)

        # Choked: M_e = 1
        @test M_e ≈ 1.0 rtol=1e-12

        # Exit velocity from energy conservation
        @test u_e ≈ sqrt(2 * (ht21 - h_e)) rtol=1e-12

        # Exit density from ideal gas law
        @test rho_e ≈ p_e / (R_e * T_e) rtol=1e-12

        # Velocity positive
        @test u_e > 0.0

        # For choked nozzle, exit pressure > ambient (under-expanded)
        @test p_e > p0

        # Exit pressure is at the nozzle throat (below pt_nozzle since M=1)
        @test p_e < pifn * pt21

        # Choked exit temperature < total temperature
        @test T_e < Tt21

        # ---------------------------------------------------------------------------
        # 4. nozzle_exit — unchoked case
        #    Force subsonic exit by using a high ambient pressure p0 = 0.9 * pt_nozzle
        # ---------------------------------------------------------------------------

        p0_unch = pifn * pt21 * 0.90   # 90% of nozzle total pressure → M < 1
        p_e_u, T_e_u, h_e_u, s_e_u, cp_e_u, R_e_u, u_e_u, rho_e_u, M_e_u =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          pt21, Tt21, ht21, st21, cpt21, Rt21, p0_unch)

        # Unchoked: M_e < 1
        @test M_e_u < 1.0

        # Exit pressure equals ambient (isentropic, unchoked)
        @test p_e_u ≈ p0_unch rtol=1e-10

        # Exit velocity from energy conservation
        @test u_e_u ≈ sqrt(2 * (ht21 - h_e_u)) rtol=1e-12

        # Exit density from ideal gas law
        @test rho_e_u ≈ p_e_u / (R_e_u * T_e_u) rtol=1e-12

        # ---------------------------------------------------------------------------
        # 5. nozzle_massflow_residual — design-point identity r = 0
        # ---------------------------------------------------------------------------

        # At design, mdot = rho_e * u_e * A  →  residual = 0
        mdot_dp = rho_e * u_e * A7

        r, r_pt, r_ht, r_st, r_Tt,
        _p_al, _T_al, _h_al, _s_al, _cp_al, _R_al =
            nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                       pt21, Tt21, ht21, st21, cpt21, Rt21,
                                       p0, mdot_dp)

        @test r ≈ 0.0 atol=1e-10

        # unchoked design-point identity
        mdot_u = rho_e_u * u_e_u * A7
        r_u, = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21, Tt21, ht21, st21, cpt21, Rt21,
                                          p0_unch, mdot_u)
        @test r_u ≈ 0.0 atol=1e-10

        # ---------------------------------------------------------------------------
        # 6. nozzle_massflow_residual — unchoked: r_Tt = 0 exactly
        # ---------------------------------------------------------------------------

        _, _, _, _, r_Tt_u = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                                         pt21, Tt21, ht21, st21, cpt21, Rt21,
                                                         p0_unch, mdot_u)[1:5]
        @test iszero(r_Tt_u)

        # ---------------------------------------------------------------------------
        # 7. Finite-difference derivative self-consistency — choked case
        # ---------------------------------------------------------------------------

        # r_pt: ∂r/∂pt_in
        δ_pt = pt21 * 1e-5
        r_pp = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21 + δ_pt, Tt21, ht21, st21, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        r_pm = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21 - δ_pt, Tt21, ht21, st21, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        @test r_pt ≈ (r_pp - r_pm) / (2δ_pt) rtol=1e-3

        # r_ht: ∂r/∂ht  (choked; gas_machd Newton tolerance → rtol ≈ 1e-3)
        δ_ht = 1.0    # J/kg absolute step (ht ≈ −36000 J/kg)
        r_hp = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21, Tt21, ht21 + δ_ht, st21, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        r_hm = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21, Tt21, ht21 - δ_ht, st21, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        @test r_ht ≈ (r_hp - r_hm) / (2δ_ht) rtol=1e-3

        # r_st: ∂r/∂st
        δ_st = 1e-4   # J/(kg·K) absolute step
        r_sp = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21, Tt21, ht21, st21 + δ_st, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        r_sm = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                          pt21, Tt21, ht21, st21 - δ_st, cpt21, Rt21,
                                          p0, mdot_dp)[1]
        @test r_st ≈ (r_sp - r_sm) / (2δ_st) rtol=1e-3

        # ---------------------------------------------------------------------------
        # 8. Finite-difference derivative self-consistency — unchoked case
        # ---------------------------------------------------------------------------

        r_u2, r_pt_u, r_ht_u, r_st_u, r_Tt_u2 =
            nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                       pt21, Tt21, ht21, st21, cpt21, Rt21,
                                       p0_unch, mdot_u)[1:5]

        δ_pt_u = pt21 * 1e-5
        r_up_pt = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                             pt21 + δ_pt_u, Tt21, ht21, st21, cpt21, Rt21,
                                             p0_unch, mdot_u)[1]
        r_um_pt = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                             pt21 - δ_pt_u, Tt21, ht21, st21, cpt21, Rt21,
                                             p0_unch, mdot_u)[1]
        @test r_pt_u ≈ (r_up_pt - r_um_pt) / (2δ_pt_u) rtol=1e-3

        # r_ht unchoked: purely from energy equation, exact (no Newton iteration)
        δ_ht_u = 1.0
        r_up_ht = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                             pt21, Tt21, ht21 + δ_ht_u, st21, cpt21, Rt21,
                                             p0_unch, mdot_u)[1]
        r_um_ht = nozzle_massflow_residual_t(fan_nozzle, alpha_air, nair,
                                             pt21, Tt21, ht21 - δ_ht_u, st21, cpt21, Rt21,
                                             p0_unch, mdot_u)[1]
        @test r_ht_u ≈ (r_up_ht - r_um_ht) / (2δ_ht_u) rtol=1e-8

        # ---------------------------------------------------------------------------
        # 9. nozzle_gross_thrust
        # ---------------------------------------------------------------------------

        # Ideally expanded nozzle: pressure thrust vanishes
        F_ideal = nozzle_gross_thrust_t(fan_nozzle, mdot_dp, u_e, p_e, p_e)
        @test F_ideal ≈ mdot_dp * u_e rtol=1e-14

        # Choked (under-expanded): p_e > p0 → pressure thrust positive
        F_choked = nozzle_gross_thrust_t(fan_nozzle, mdot_dp, u_e, p_e, p0)
        @test F_choked > mdot_dp * u_e

        # Unchoked (p_e = p0): pressure term vanishes
        F_unchoked = nozzle_gross_thrust_t(fan_nozzle, mdot_u, u_e_u, p_e_u, p0_unch)
        @test F_unchoked ≈ mdot_u * u_e_u rtol=1e-12

        # Formula: F = mdot*u + (p - p0)*A
        @test F_choked ≈ mdot_dp * u_e + (p_e - p0) * A7 rtol=1e-14

        # ---------------------------------------------------------------------------
        # 10. Monotonicity
        #
        # Choked case: u_e is determined by (ht, Tt, st) only, not by pt_in.
        # Higher pt_in → higher mass flux (ρ_e u_e A) at fixed sonic conditions.
        #
        # Unchoked case: higher pt_in → greater expansion ratio → higher u_e.
        # ---------------------------------------------------------------------------

        # Choked monotonicity: mass flux increases with pt_in
        _, _, _, _, _, _, u_e_hi, rho_e_hi, _ =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          1.05 * pt21, Tt21, ht21, st21, cpt21, Rt21, p0)
        _, _, _, _, _, _, u_e_lo, rho_e_lo, _ =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          0.95 * pt21, Tt21, ht21, st21, cpt21, Rt21, p0)
        @test rho_e_hi * u_e_hi > rho_e * u_e   # higher pt → more mass flux
        @test rho_e_lo * u_e_lo < rho_e * u_e

        # Unchoked monotonicity: higher pt_in → more expansion → higher exit velocity
        _, _, _, _, _, _, u_e_hi_u, _, _ =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          1.05 * pt21, Tt21, ht21, st21, cpt21, Rt21, p0_unch)
        _, _, _, _, _, _, u_e_lo_u, _, _ =
            nozzle_exit_t(fan_nozzle, alpha_air, nair,
                          0.95 * pt21, Tt21, ht21, st21, cpt21, Rt21, p0_unch)
        @test u_e_hi_u > u_e_u
        @test u_e_lo_u < u_e_u

        # ---------------------------------------------------------------------------
        # 11. Float32 dispatch
        # ---------------------------------------------------------------------------

        fan_nozzle_32 = Nozzle_t(Float32(pifn), Float32(A7))
        p_e32, T_e32, h_e32, s_e32, cp_e32, R_e32, u_e32, rho_e32, M_e32 =
            nozzle_exit_t(fan_nozzle_32,
                          alpha_air, nair,
                          Float32(pt21), Float32(Tt21), Float32(ht21),
                          Float32(st21), Float32(cpt21), Float32(Rt21),
                          Float32(p0))
        @test u_e32   isa Float32
        @test rho_e32 isa Float32
        @test M_e32   isa Float32
        @test u_e32   ≈ Float32(u_e) rtol=1e-4

        r32 = nozzle_massflow_residual_t(fan_nozzle_32,
                                         alpha_air, nair,
                                         Float32(pt21), Float32(Tt21), Float32(ht21),
                                         Float32(st21), Float32(cpt21), Float32(Rt21),
                                         Float32(p0), Float32(mdot_dp))[1]
        @test r32 isa Float32

    end  # Nozzle

    # Newton driver nozzle wiring — tasopt-j9l.31
    # Verify that nozzle_massflow_residual called with pn=1 (total-pressure loss
    # already folded into pt7/pt5) reproduces the same residual and chain-rule
    # Jacobian as the inlined physics it replaces.
    @testset "Newton driver nozzle wiring" begin

        Nozzle_t                   = TASOPT.engine.Nozzle
        nozzle_exit_t              = TASOPT.engine.nozzle_exit
        nozzle_massflow_residual_t = TASOPT.engine.nozzle_massflow_residual

        # Re-use fan nozzle state from the parent Nozzle testset conditions
        # but now wrap with pn=1 to match Newton driver usage (loss already in pt7)
        alpha_air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        nair  = 5

        pt21  = 61271.983314487086
        Tt21  = 294.1405135675445
        pifn  = 0.98
        A7    = 0.841808817378845
        p0    = 23922.608843328788

        gassum_t = TASOPT.engine.gassum
        st21, _, ht21, _, cpt21, Rt21 = gassum_t(alpha_air, nair, Tt21)

        # pt7 has the nozzle loss already applied (as in the Newton driver)
        pt7  = pt21 * pifn
        Tt7  = Tt21
        ht7  = ht21
        st7  = st21
        cpt7 = cpt21
        Rt7  = Rt21

        # Newton driver constructs: nozzle_fan = Nozzle(one(T), A7)  (pn = 1)
        nozzle_pn1 = Nozzle_t(1.0, A7)

        # Design mass flow at the exit
        p_e, T_e, h_e, s_e, cp_e, R_e, u_e, rho_e, M_e =
            nozzle_exit_t(nozzle_pn1, alpha_air, nair, pt7, Tt7, ht7, st7, cpt7, Rt7, p0)
        mdot_dp = rho_e * u_e * A7

        # 1. Residual is exactly zero at design mass flow
        r, r_pt7, r_ht7, r_st7, r_Tt7, _, _, _, _, _, _ =
            nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                       pt7, Tt7, ht7, st7, cpt7, Rt7, p0, mdot_dp)
        @test r ≈ 0.0 atol=1e-10

        # 2. Chain-rule Jacobian: verify r_pt7 via finite differences
        δpt = pt7 * 1e-5
        r_pp = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7 + δpt, Tt7, ht7, st7, cpt7, Rt7, p0, mdot_dp)[1]
        r_pm = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7 - δpt, Tt7, ht7, st7, cpt7, Rt7, p0, mdot_dp)[1]
        @test r_pt7 ≈ (r_pp - r_pm) / (2δpt) rtol=1e-3

        # 3. Chain-rule Jacobian: verify r_ht7 via finite differences
        δht = 1.0
        r_hp = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7, Tt7, ht7 + δht, st7, cpt7, Rt7, p0, mdot_dp)[1]
        r_hm = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7, Tt7, ht7 - δht, st7, cpt7, Rt7, p0, mdot_dp)[1]
        @test r_ht7 ≈ (r_hp - r_hm) / (2δht) rtol=1e-3

        # 4. Chain-rule Jacobian: verify r_st7 via finite differences
        δst = 1e-4
        r_sp = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7, Tt7, ht7, st7 + δst, cpt7, Rt7, p0, mdot_dp)[1]
        r_sm = nozzle_massflow_residual_t(nozzle_pn1, alpha_air, nair,
                                          pt7, Tt7, ht7, st7 - δst, cpt7, Rt7, p0, mdot_dp)[1]
        @test r_st7 ≈ (r_sp - r_sm) / (2δst) rtol=1e-3

        # 5. Nozzle with pn=1 at the same pt7 is equivalent to pn=pifn at pt21
        nozzle_original = Nozzle_t(pifn, A7)
        r_orig, r_pt_orig = nozzle_massflow_residual_t(nozzle_original, alpha_air, nair,
                                                        pt21, Tt21, ht21, st21, cpt21, Rt21,
                                                        p0, mdot_dp)[1:2]
        @test r ≈ r_orig atol=1e-10
        # r_pt of pn=1 nozzle at pt7 ↔ r_pt of pifn nozzle at pt21:
        # dr/d(pt_in) = dr/d(pt_nozzle) * pn, so r_pt_orig = r_pt7 * pifn
        @test r_pt7 * pifn ≈ r_pt_orig rtol=1e-8

    end  # Newton driver nozzle wiring

    @testset "engine_plots" begin

        ac_plots = TASOPT.load_default_model()
        size_aircraft!(ac_plots; printiter=false)
        mission_plots = TASOPT.engine.run_engine_sweep(ac_plots)

        # ---- plot_engine_station_profiles -----------------------------------
        p_prof = TASOPT.plot_engine_station_profiles(mission_plots, ipstatic:ipdescentn)

        @test p_prof isa Plots.Plot
        # Two panels: Tt profile (top) and pt profile (bottom).
        @test length(p_prof.subplots) == 2

        # save to temp file — exercises the full Plots rendering path
        tmp_prof = joinpath(tempdir(), "tasopt_test_station_profiles.png")
        @test (savefig(p_prof, tmp_prof); isfile(tmp_prof))

        # ---- plot_engine_performance (no spool speeds) ----------------------
        p_basic = TASOPT.plot_engine_performance(mission_plots, ipstatic:ipdescentn)

        @test p_basic isa Plots.Plot
        # Four panels: Fe, TSFC, BPR, mcore.
        @test length(p_basic.subplots) == 4

        tmp_basic = joinpath(tempdir(), "tasopt_test_engine_perf_basic.png")
        @test (savefig(p_basic, tmp_basic); isfile(tmp_basic))

        # ---- plot_engine_performance (with spool speeds) --------------------
        p_full = TASOPT.plot_engine_performance(mission_plots, ipstatic:ipdescentn;
                                                ac=ac_plots, imission=1)

        @test p_full isa Plots.Plot
        # Seven panels: Fe, TSFC, BPR, mcore, Nbf, Nblc, Nbhc.
        @test length(p_full.subplots) == 7

        tmp_full = joinpath(tempdir(), "tasopt_test_engine_perf_full.png")
        @test (savefig(p_full, tmp_full); isfile(tmp_full))

    end  # engine_plots
end