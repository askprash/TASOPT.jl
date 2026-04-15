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
        ac.parg[igWMTO] = MTOW
        ac.parg[igfeng] = feng

        TASOPT.engine.constant_TSFC_engine!(ac, 0, 1, ip, 0)

        @test ac.pare[iemfuel,ip,1] ≈ neng*TSFC*Fe/gee rtol = 1e-10

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

        pare_copy = copy(collect(pare_orig))
        TASOPT.engine.design_state_to_pare!(ds, pare_copy)

        for idx in [ieA2, ieA25, ieA5, ieA7,
                    ieNbfD, ieNblcD, ieNbhcD, ieNbhtD, ieNbltD,
                    iembfD, iemblcD, iembhcD, iembhtD, iembltD,
                    iepifD, iepilcD, iepihcD, iepihtD, iepiltD]
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
    # run_engine_sweep / SweepResult / write_sweep_csv
    # ======================================================================
    @testset "run_engine_sweep" begin
        # ------------------------------------------------------------------
        # Setup: use the same sized aircraft from the design-point test.
        # ------------------------------------------------------------------
        ac = TASOPT.load_default_model()
        size_aircraft!(ac; printiter=false)

        # ------------------------------------------------------------------
        # Sweep all regular mission points (ipstatic:ipdescentn = 1:16).
        # After sizing, every pare column is already converged, so the sweep
        # should reproduce those converged values.
        # ------------------------------------------------------------------
        sweep = TASOPT.engine.run_engine_sweep(ac)

        # ---- Structural invariants ----
        n_pts = ipdescentn - ipstatic + 1   # should be 16
        @test length(sweep.ip_indices) == n_pts
        @test length(sweep.ip_labels)  == n_pts
        @test length(sweep.engines)    == n_pts
        @test length(sweep.Fe)         == n_pts
        @test length(sweep.TSFC)       == n_pts
        @test length(sweep.BPR)        == n_pts
        @test length(sweep.mcore)      == n_pts
        @test length(sweep.mdotf)      == n_pts

        # First and last ip indices match ipstatic and ipdescentn
        @test sweep.ip_indices[1]   == ipstatic
        @test sweep.ip_indices[end] == ipdescentn

        # ---- All performance scalars are positive (physical) ----
        @test all(x -> x > 0.0, sweep.Fe)
        @test all(x -> x > 0.0, sweep.TSFC)
        @test all(x -> x > 0.0, sweep.BPR)
        @test all(x -> x > 0.0, sweep.mcore)
        @test all(x -> x > 0.0, sweep.mdotf)

        # ---- Mach and altitude are populated for all points ----
        @test all(x -> x >= 0.0, sweep.Mach)
        @test all(x -> x >= 0.0, sweep.alt)

        # ---- Engine states match pare at every mission point ----
        # cruise1 index within the sweep vector
        k_cr = ipcruise1 - ipstatic + 1
        eng_cr = sweep.engines[k_cr]
        pare_cr = view(ac.pare, :, ipcruise1, 1)

        @test eng_cr.Tt4  ≈ pare_cr[ieTt4]   rtol = 1e-12
        @test eng_cr.pt3  ≈ pare_cr[iept3]   rtol = 1e-12
        @test eng_cr.Tt49 ≈ pare_cr[ieTt49]  rtol = 1e-12
        @test sweep.BPR[k_cr] ≈ pare_cr[ieBPR]    rtol = 1e-12
        @test sweep.Fe[k_cr]  ≈ pare_cr[ieFe]     rtol = 1e-12
        @test sweep.TSFC[k_cr] ≈ pare_cr[ieTSFC]  rtol = 1e-12
        @test sweep.mcore[k_cr] ≈ pare_cr[iemcore] rtol = 1e-12
        @test sweep.mdotf[k_cr] ≈ pare_cr[iemfuel] rtol = 1e-12

        # Check a climb point too (ipclimb1 uses FixedTt4OffDes mode)
        k_cl = ipclimb1 - ipstatic + 1
        eng_cl = sweep.engines[k_cl]
        pare_cl = view(ac.pare, :, ipclimb1, 1)

        @test eng_cl.Tt4  ≈ pare_cl[ieTt4]   rtol = 1e-12
        @test eng_cl.pt3  ≈ pare_cl[iept3]   rtol = 1e-12

        # ---- Thermodynamic invariants on cruise point ----
        # Total temperature must rise through compressor, fall through turbine
        @test eng_cr.Tt25 > eng_cr.Tt19   # LPC adds work
        @test eng_cr.Tt3  > eng_cr.Tt25   # HPC adds work
        @test eng_cr.Tt4  > eng_cr.Tt3    # combustor
        @test eng_cr.Tt41 < eng_cr.Tt4    # turbine inlet (cooling dilution)
        @test eng_cr.Tt45 < eng_cr.Tt41   # HPT
        @test eng_cr.Tt49 < eng_cr.Tt45   # LPT

        # ---- Custom ip_range: only the two cruise points ----
        sweep_cr2 = TASOPT.engine.run_engine_sweep(ac;
                        ip_range = ipcruise1:ipcruise2)
        @test length(sweep_cr2.ip_indices) == 2
        @test sweep_cr2.ip_indices[1] == ipcruise1
        @test sweep_cr2.ip_indices[2] == ipcruise2
        @test all(x -> x > 0.0, sweep_cr2.Fe)

        # ---- CSV serialisation round-trip ----
        buf = IOBuffer()
        TASOPT.engine.write_sweep_csv(buf, sweep_cr2)
        csv_str = String(take!(buf))

        # Header row must contain expected column names
        header_line = split(csv_str, "\n")[1]
        @test occursin("ip",        header_line)
        @test occursin("label",     header_line)
        @test occursin("Fe_N",      header_line)
        @test occursin("TSFC_kg_Ns",  header_line)
        @test occursin("BPR",         header_line)
        @test occursin("mfuel_kg_s",  header_line)
        @test occursin("Tt4_K",     header_line)
        @test occursin("pt3_Pa",    header_line)

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
        sweep_cr2 = TASOPT.engine.run_engine_sweep(ac; ip_range=ipcruise1:ipcruise2)

        # ------------------------------------------------------------------
        # write_sweep_toml writes valid, parseable TOML.
        # ------------------------------------------------------------------
        buf = IOBuffer()
        TASOPT.engine.write_sweep_toml(buf, sweep_cr2)
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

        # Stations not populated from pare have zero total temperature
        @test stations_p1["st19c"]["Tt"] == 0.0
        @test stations_p1["st25c"]["Tt"] == 0.0
        @test stations_p1["st4a"]["Tt"]  == 0.0
        @test stations_p1["st49c"]["Tt"] == 0.0

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
        sweep = TASOPT.engine.run_engine_sweep(ac)

        baseline_path = TASOPT.engine.ENGINE_BASELINE_PATH
        @test isfile(baseline_path)
        baseline = TOML.parsefile(baseline_path)

        n_pts = length(sweep.ip_indices)
        baseline_pts = baseline["points"]
        @test length(baseline_pts) == n_pts

        # ---- Comparison tolerance: one part per trillion ----
        rtol_tol = 1e-12

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
            eng = sweep.engines[k]

            @test sweep.ip_labels[k] == lbl  # ordering sanity check

            # Top-level performance scalars
            chk!(sweep.Fe[k],    bp["Fe_N"],       "[$lbl] Fe_N")
            chk!(sweep.TSFC[k],  bp["TSFC_kg_Ns"], "[$lbl] TSFC_kg_Ns")
            chk!(sweep.BPR[k],   bp["BPR"],        "[$lbl] BPR")
            chk!(sweep.mcore[k], bp["mcore_kg_s"], "[$lbl] mcore_kg_s")
            chk!(sweep.mdotf[k], bp["mfuel_kg_s"], "[$lbl] mfuel_kg_s")
            chk!(sweep.alt[k],   bp["alt_m"],      "[$lbl] alt_m")
            chk!(sweep.Mach[k],  bp["Mach"],       "[$lbl] Mach")
            chk!(eng.M0,         bp["M0"],         "[$lbl] M0")
            chk!(eng.T0,         bp["T0_K"],       "[$lbl] T0_K")
            chk!(eng.p0,         bp["p0_Pa"],      "[$lbl] p0_Pa")
            chk!(eng.a0,         bp["a0_m_s"],     "[$lbl] a0_m_s")

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
end