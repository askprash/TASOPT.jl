# Equivalence of the IdealGasThermo.jl-based gas properties in `TASOPT.engine`
# with the table-based implementation they replaced.
#
# The reference is a verbatim copy of the old code in test/legacy_thermo/.
# The two implementations use different datums for h and s (IdealGasThermo
# uses formation-inclusive h and absolute standard entropy, the legacy tables
# used their own offsets), so h and s are compared as changes from 298.15 K.
# Only differences in h and s enter the cycle calculations at fixed
# composition.

isdefined(@__MODULE__, :LegacyThermo) ||
    include(joinpath(@__DIR__, "legacy_thermo", "LegacyThermo.jl"))

@testset "Thermo equivalence with legacy tables" verbose = true begin
    E = TASOPT.engine
    L = LegacyThermo
    Tref = 298.15
    Tgrid = 200.0:10.0:2500.0

    # Largest relative deviations of new vs legacy over Tgrid:
    #   cp: |cp_new/cp_old - 1|
    #   dh: |Δh_new - Δh_old| / (cp_old*T),  Δh = h(T) - h(298.15)
    #   ds: |Δs_new - Δs_old| / R_old,       Δs = s(T) - s(298.15)
    #   R : |R_new/R_old - 1|
    function max_deviation(props_new, props_old)
        s_n0, _, h_n0, _, _, _ = props_new(Tref)
        s_o0, _, h_o0, _, _, _ = props_old(Tref)
        dcp = ddh = dds = dR = 0.0
        for T in Tgrid
            s_n, _, h_n, _, cp_n, r_n = props_new(T)
            s_o, _, h_o, _, cp_o, r_o = props_old(T)
            dcp = max(dcp, abs(cp_n / cp_o - 1))
            ddh = max(ddh, abs((h_n - h_n0) - (h_o - h_o0)) / (cp_o * T))
            dds = max(dds, abs((s_n - s_n0) - (s_o - s_o0)) / r_o)
            dR = max(dR, abs(r_n / r_o - 1))
        end
        return (cp = dcp, dh = ddh, ds = dds, R = dR)
    end

    @testset "Pure species, 200-2500 K" begin
        # (igas, name, cp, dh, ds, R) tolerances.
        # N2, O2, CO2, Ar: tables and NASA-9 fits agree to ~0.1%.
        # H2O: legacy table cp falls below the NASA data above ~1500 K
        #      (-1.8% at 2500 K); the NASA-9 fit is the better source.
        # H2:  legacy table cp is up to 0.6% low near 200 K.
        tolerances = [
            (1, "N2", 2e-3, 1e-3, 5e-3, 1e-3),
            (2, "O2", 2e-3, 1e-3, 5e-3, 1e-3),
            (3, "CO2", 2e-3, 1e-3, 5e-3, 1e-3),
            (5, "Ar", 2e-3, 1e-3, 5e-3, 1e-3),
            (4, "H2O", 2e-2, 1e-2, 6e-2, 1e-3),
            (40, "H2", 1e-2, 1e-3, 5e-2, 1e-3),
        ]
        for (igas, name, tcp, tdh, tds, tR) in tolerances
            dev = max_deviation(T -> E.gasfun(igas, T), T -> L.gasfun(igas, T))
            @testset "$name" begin
                @test dev.cp < tcp
                @test dev.dh < tdh
                @test dev.ds < tds
                @test dev.R < tR
            end
        end
    end

    @testset "Mixtures, 200-2500 K" begin
        air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        n = length(air)
        products(far) = (air .+ far .* E.gasfuel(24, n + 1)[1:n]) ./ (1 + far)
        mixtures = [
            ("air", air, 1e-3, 2e-4, 1e-3, 5e-4),
            ("products FAR=0.01", products(0.01), 2.5e-3, 1e-3, 5e-3, 5e-4),
            ("products FAR=0.02", products(0.02), 2.5e-3, 1e-3, 5e-3, 5e-4),
            ("products FAR=0.03", products(0.03), 2.5e-3, 1e-3, 5e-3, 5e-4),
        ]
        for (name, alpha, tcp, tdh, tds, tR) in mixtures
            dev = max_deviation(T -> E.gassum(alpha, n, T), T -> L.gassum(alpha, n, T))
            @testset "$name" begin
                @test dev.cp < tcp
                @test dev.dh < tdh
                @test dev.ds < tds
                @test dev.R < tR
            end
        end
    end

    @testset "Fuels" begin
        # The legacy fuel tables are coarse (cp off by up to 35%), so fuels
        # are not compared point-wise. What the cycle uses is the heat
        # released by burning them, so that is what must be preserved.
        fuels = [(11, "CH4"), (12, "C2H6"), (13, "C3H8"), (14, "C4H10"),
                 (18, "C8H18"), (24, "C14H30"), (40, "H2")]
        wchon = [12.01078, 1.00795, 15.99943, 14.00672] # as in gasfuel
        for (igas, name) in fuels
            @testset "$name" begin
                @test E.fuelLHV(igas) ≈ L.fuelLHV(igas) rtol = 1e-3
                # thermo data and stoichiometry describe the same molecule
                MW = sum(wchon .* E.gaschem(igas))
                @test E.gasfun(igas, 300.0)[6] ≈ 8314.462618 / MW rtol = 1e-4
            end
        end
        @testset "Jet-A (C14H30) burner fuel/air ratio" begin
            air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]
            beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
            for (Tair, Tfuel, Tb) in ((800.0, 300.0, 1200.0), (850.0, 300.0, 1600.0), (900.0, 400.0, 1900.0))
                f_new, _ = E.gas_burn(air, beta, E.gasfuel(24, 6), 6, 24, Tair, Tfuel, Tb, 0.0)
                f_old, _ = L.gas_burn(air, beta, L.gasfuel(24, 6), 6, 24, Tair, Tfuel, Tb, 0.0)
                @test f_new ≈ f_old rtol = 1e-3
            end
        end
    end

    @testset "H2 below 200 K (LH2 coolant range)" begin
        # NASA-9 fits start at 200 K. Below that, H2 must follow the legacy
        # tabulated data (20.369-200 K), which the heat-exchanger models rely
        # on for LH2 coolant and the burner for cold H2 fuel.
        igas = 40
        Tlow = 200.0
        _, _, h_n200, _, _, _ = E.gasfun(igas, Tlow)
        s_n200 = E.gasfun(igas, Tlow)[1]
        s_o200, _, h_o200, _, _, _ = L.gasfun(igas, Tlow)
        for T in (20.369, 25.0, 40.0, 60.0, 100.0, 150.0, 190.0)
            s_n, s_tn, h_n, h_tn, cp_n, r_n = E.gasfun(igas, T)
            s_o, _, h_o, _, cp_o, _ = L.gasfun(igas, T)
            @test cp_n ≈ cp_o rtol = 1e-3
            @test abs((h_n - h_n200) - (h_o - h_o200)) < 1e-3 * cp_o * T
            @test abs((s_n - s_n200) - (s_o - s_o200)) < 1e-3 * r_n
            @test h_tn == cp_n
            @test s_tn ≈ cp_n / T rtol = 1e-14
        end
        # h and s are continuous where the data source switches
        ε = 1e-9
        _, _, h_lo, _, cp_lo, _ = E.gasfun(igas, Tlow - ε)
        _, _, h_hi, _, cp_hi, _ = E.gasfun(igas, Tlow + ε)
        @test abs(h_hi - h_lo) < 1e-6 * cp_hi
        @test abs(E.gasfun(igas, Tlow + ε)[1] - E.gasfun(igas, Tlow - ε)[1]) < 1e-6
        # gasfund (summed by gassum/gassumd) agrees with gasfun below 200 K,
        # and its dcp/dT matches a central difference of cp
        for T in (30.0, 120.0)
            s, s_t, h, h_t, cp, cp_t, r = E.gasfund(igas, T)
            @test (s, s_t, h, h_t, cp, r) == E.gasfun(igas, T)
            @test cp_t ≈ (E.gasfun(igas, T + 0.01)[5] - E.gasfun(igas, T - 0.01)[5]) / 0.02 rtol = 1e-3
        end
        # enthalpy inversion works in the cryogenic range
        for T in (25.0, 60.0, 150.0)
            hT = E.gasfun(igas, T)[3]
            @test E.gas_tset_single(igas, hT, 100.0) ≈ T rtol = 1e-8
        end
    end

    @testset "Process functions agree with legacy" begin
        air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        n = length(air)
        lam = (air .+ 0.025 .* E.gasfuel(24, n + 1)[1:n]) ./ 1.025 # burned gas
        state(M, a, T) = (s = M.gassum(a, n, T); (s[1], s[3], s[5], s[6])) # s, h, cp, R
        rel(a, b) = abs(a / b - 1)

        @testset "gas_tset" begin
            # temperature reached after adding Δh to a gas at T1
            for (a, T1, dh) in ((air, 250.0, 5e4), (air, 600.0, 3e5), (lam, 1500.0, -3e5))
                T_new = E.gas_tset(a, n, E.gassum(a, n, T1)[3] + dh, T1)
                T_old = L.gas_tset(a, n, L.gassum(a, n, T1)[3] + dh, T1)
                @test rel(T_new, T_old) < 5e-4
            end
        end
        @testset "gas_prat" begin
            for (a, T1, pr, ep) in ((air, 288.15, 1.6, 0.9), (air, 400.0, 10.0, 0.9),
                                    (lam, 1600.0, 0.2, 1 / 0.9), (lam, 1100.0, 0.5, 1 / 0.9))
                sn, hn, cpn, rn = state(E, a, T1)
                so, ho, cpo, ro = state(L, a, T1)
                p_n, T_n, h_n = E.gas_prat(a, n, 1e5, T1, hn, sn, cpn, rn, pr, ep)
                p_o, T_o, h_o = L.gas_prat(a, n, 1e5, T1, ho, so, cpo, ro, pr, ep)
                @test p_n == p_o
                @test rel(T_n, T_o) < 5e-4
                @test rel(h_n - hn, h_o - ho) < 1e-3 # work
            end
        end
        @testset "gas_delh" begin
            for (a, T1, dh, ep) in ((air, 288.15, 3e5, 0.9), (air, 500.0, 5e5, 0.9), (lam, 1600.0, -4e5, 1 / 0.9))
                sn, hn, cpn, rn = state(E, a, T1)
                so, ho, cpo, ro = state(L, a, T1)
                p_n, T_n = E.gas_delh(a, n, 1e5, T1, hn, sn, cpn, rn, dh, ep)
                p_o, T_o = L.gas_delh(a, n, 1e5, T1, ho, so, cpo, ro, dh, ep)
                @test rel(T_n, T_o) < 5e-4
                @test rel(p_n, p_o) < 2e-3
            end
        end
        @testset "gas_mach" begin
            for (a, T1, M) in ((air, 300.0, 0.5), (air, 300.0, 0.9), (lam, 1200.0, 1.0))
                sn, hn, cpn, rn = state(E, a, T1)
                so, ho, cpo, ro = state(L, a, T1)
                p_n, T_n = E.gas_mach(a, n, 1e5, T1, hn, sn, cpn, rn, 0.0, M, 1.0)
                p_o, T_o = L.gas_mach(a, n, 1e5, T1, ho, so, cpo, ro, 0.0, M, 1.0)
                @test rel(T_n, T_o) < 5e-4
                @test rel(p_n, p_o) < 5e-4
            end
        end
        @testset "gas_mass" begin
            for (a, T1, mflux) in ((air, 300.0, 150.0), (lam, 1200.0, 60.0))
                sn, hn, cpn, rn = state(E, a, T1)
                so, ho, cpo, ro = state(L, a, T1)
                p_n, T_n = E.gas_mass(a, n, 1e5, T1, hn, sn, cpn, rn, mflux, 0.3)
                p_o, T_o = L.gas_mass(a, n, 1e5, T1, ho, so, cpo, ro, mflux, 0.3)
                @test rel(T_n, T_o) < 5e-4
                @test rel(p_n, p_o) < 5e-4
            end
        end
        @testset "gas_burn" begin
            # H2 makes 8.9 kg of H2O per kg of fuel, so its fuel/air ratio
            # carries the H2O sensible-enthalpy difference (legacy table vs
            # NASA-9: +0.29% from 298 to 1600 K): f +0.15%. Jet-A +0.07%,
            # CH4 +0.02%.
            air6 = [air; 0.0]
            beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
            for (ifuel, Tair, Tfuel, Tb) in ((24, 800.0, 300.0, 1600.0), (40, 800.0, 300.0, 1600.0),
                                             (40, 800.0, 50.0, 1600.0), (11, 700.0, 250.0, 1400.0))
                f_n, l_n = E.gas_burn(air6, beta, E.gasfuel(ifuel, 6), 6, ifuel, Tair, Tfuel, Tb, 0.0)
                f_o, l_o = L.gas_burn(air6, beta, L.gasfuel(ifuel, 6), 6, ifuel, Tair, Tfuel, Tb, 0.0)
                @test rel(f_n, f_o) < 2e-3
                @test maximum(abs.(l_n .- l_o)) < 2e-4
            end
        end
    end

    @testset "Partials helpers match gas_pratd / gas_delhd" begin
        # tfoper uses the NamedTuple partials functions and per-direction
        # composition sensitivities; they must reproduce the array-based
        # gas_pratd / gas_delhd exactly (same property data on both sides).
        n = 5
        lam = [0.74, 0.15, 0.06, 0.03, 0.02]
        dirs = ([0.01, -0.03, 0.02, 0.005, -0.005], [0.0, 0.001, -0.002, 0.003, -0.002],
                [-0.02, 0.01, 0.0, 0.01, 0.0], [0.1, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 1.0], [0.3, -0.1, 0.2, -0.1, 0.0], [0.0, 1.0, 0.0, 0.0, 0.0])
        s1, _, h1, _, cp1, r1 = E.gassum(lam, n, 1500.0)
        args = (lam, n, 1.0e6, 1500.0, h1, s1, cp1, r1)
        # relative 1e-10, with an absolute floor of 1e-10 for quantities that
        # are exactly zero in one form and cancel to round-off in the other
        # (e.g. dh/dalpha at fixed delh)
        close(a, b) = abs(a - b) <= 1e-10 * max(abs(a), abs(b), 1.0)

        o = E.gas_pratd(args..., 0.3, 1 / 0.9)
        st, d = E.gas_pressure_ratio_partials(args..., 0.3, 1 / 0.9)
        @test all(close.((st.p, st.T, st.h, st.s, st.cp, st.R), o[1:6]))
        @test all(close.((d.p_po, d.p_so, d.T_so, d.h_so, d.s_so, d.p_pratio, d.T_pratio, d.h_pratio,
                          d.s_pratio, d.p_epol, d.T_epol, d.h_epol, d.s_epol), o[7:19]))
        p_al, T_al, h_al, s_al, cp_al, R_al = o[20:25]
        for (k, dk) in enumerate(E.gas_pressure_ratio_composition_directions(lam, n, st.T, s1, dirs...))
            @test close(dk.T, sum(T_al .* dirs[k]))
            @test close(dk.h, sum(h_al .* dirs[k]))
            @test close(dk.s, sum(s_al .* dirs[k]))
            @test close(dk.cp, sum(cp_al .* dirs[k]))
            @test close(dk.R, sum(R_al .* dirs[k]))
        end

        o = E.gas_delhd(args..., -3e5, 1 / 0.9)
        st, d = E.gas_enthalpy_change_partials(args..., -3e5, 1 / 0.9)
        @test all(close.((st.p, st.T, st.h, st.s, st.cp, st.R), o[1:6]))
        @test all(close.((d.p_so, d.p_po, d.p_epol, d.p_ho, d.T_ho, d.h_ho, d.s_ho,
                          d.p_delh, d.T_delh, d.h_delh, d.s_delh), o[7:17]))
        p_al, T_al, h_al, s_al, cp_al, R_al = o[18:23]
        for (k, dk) in enumerate(E.gas_enthalpy_change_composition_directions(lam, n, st.p, st.T, s1, 1 / 0.9, dirs...))
            @test close(dk.p, sum(p_al .* dirs[k]))
            @test close(dk.T, sum(T_al .* dirs[k]))
            @test close(dk.h, sum(h_al .* dirs[k]))
            @test close(dk.s, sum(s_al .* dirs[k]))
        end
    end

    @testset "Precomputed mixtures match the species sums" begin
        # gassum/gassumd accept a GasMixture in place of the mass-fraction
        # vector: the NASA-9 coefficients of the constituents summed once, with
        # the mass fractions as weights. Properties are linear in the
        # coefficients, so the result must equal the species loop to round-off,
        # for any coefficient vector (including the composition directions
        # tfoper! passes, which have negative entries and do not sum to 1).
        air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        n = length(air)
        products(far) = (air .+ far .* E.gasfuel(24, n + 1)[1:n]) ./ (1 + far)
        compositions = [("air", air), ("products FAR=0.01", products(0.01)),
                        ("products FAR=0.03", products(0.03)),
                        ("direction", [0.012, -0.031, 0.02, 0.004, -0.005]),
                        ("air, 3 constituents", air[1:3])]
        Ts = [200.0, 288.15, 450.0, 999.0, 1000.0, 1001.0, 1500.0, 2500.0]
        close(a, b) = abs(a - b) <= 1e-13 * max(abs(a), abs(b), 1.0)

        for (name, alpha) in compositions
            m = length(alpha)
            mix = E.GasMixture(alpha, m)
            @testset "$name" begin
                for T in Ts
                    @test all(close.(E.gassum(mix, m, T), E.gassum(alpha, m, T)))
                    @test all(close.(E.gassumd(mix, m, T), E.gassumd(alpha, m, T)))
                end
            end
        end

        @testset "absolute entropy is comparable across compositions" begin
            # A mixture built as an IdealGasThermo FrozenGas from mole fractions
            # includes the ideal entropy of mixing, -R sum(X ln X), in s0: about
            # 169 J/kg/K for air and a different amount for products. tfoper!
            # builds composition derivatives of s from the constituent values
            # (ds/dalpha_i = s_i), which hold only for the plain mass-weighted
            # sum, so that constant must not appear.
            pmix = E.GasMixture(products(0.03), n)
            amix = E.GasMixture(air, n)
            for T in Ts
                ds_mix = E.gassum(pmix, n, T)[1] - E.gassum(amix, n, T)[1]
                ds_sum = E.gassum(products(0.03), n, T)[1] - E.gassum(air, n, T)[1]
                @test abs(ds_mix - ds_sum) < 1e-9
                @test abs(E.gassum(amix, n, T)[1] - E.gassum(air, n, T)[1]) < 1e-9
            end
        end

        @testset "process functions give the same result" begin
            lam = products(0.025)
            for (a, T1) in ((air, 288.15), (lam, 1500.0))
                mix = E.GasMixture(a, n)
                s1, _, h1, _, cp1, r1 = E.gassum(a, n, T1)
                st = (1.0e5, T1, h1, s1, cp1, r1)
                # 1e-12: the iterations amplify the 1e-16 differences of the
                # property sums (gas_mass: h is a difference of terms of order
                # cp*T, 1.7e6 J/kg for the products, and comes out 2e-13 apart)
                same(x, y) = all(abs.(x .- y) .<= 1e-12 .* max.(abs.(x), abs.(y), 1.0))
                @test E.gas_tset(mix, n, h1 + 2e5, T1) ≈ E.gas_tset(a, n, h1 + 2e5, T1) rtol = 1e-13
                @test same(E.gas_tsetd(mix, n, h1 + 2e5, T1)[1:2], E.gas_tsetd(a, n, h1 + 2e5, T1)[1:2])
                @test same(E.gas_tsetd(mix, n, h1 + 2e5, T1)[3], E.gas_tsetd(a, n, h1 + 2e5, T1)[3])
                @test same(E.gas_prat(mix, n, st..., 1.7, 0.9), E.gas_prat(a, n, st..., 1.7, 0.9))
                @test same(E.gas_delh(mix, n, st..., -2e5, 1 / 0.9), E.gas_delh(a, n, st..., -2e5, 1 / 0.9))
                @test same(E.gas_mach(mix, n, st..., 0.0, 0.8, 1.0), E.gas_mach(a, n, st..., 0.0, 0.8, 1.0))
                @test same(E.gas_mass(mix, n, st..., 100.0, 0.3), E.gas_mass(a, n, st..., 100.0, 0.3))
                pm, dm = E.gas_pressure_ratio_partials(mix, n, st..., 1.7, 0.9)
                pa, da = E.gas_pressure_ratio_partials(a, n, st..., 1.7, 0.9)
                @test same(values(pm), values(pa)) && same(values(dm), values(da))
                pm, dm = E.gas_enthalpy_change_partials(mix, n, st..., -2e5, 1 / 0.9)
                pa, da = E.gas_enthalpy_change_partials(a, n, st..., -2e5, 1 / 0.9)
                @test same(values(pm), values(pa)) && same(values(dm), values(da))
                om = E.gas_machd(mix, n, st..., 0.0, 0.8, 1.0)
                oa = E.gas_machd(a, n, st..., 0.0, 0.8, 1.0)
                @test same(om[1:21], oa[1:21])
                @test all(same(om[k], oa[k]) for k in 22:27)
                d = [0.012, -0.031, 0.02, 0.004, -0.005]
                @test same(values(E.gas_pressure_ratio_composition_direction(mix, n, 1100.0, s1, d)),
                           values(E.gas_pressure_ratio_composition_direction(a, n, 1100.0, s1, d)))
                @test same(values(E.gas_enthalpy_change_composition_direction(mix, n, 2e5, 1100.0, s1, 0.9, d)),
                           values(E.gas_enthalpy_change_composition_direction(a, n, 2e5, 1100.0, s1, 0.9, d)))
            end
        end

        @testset "no allocation, n checked, Dual numbers" begin
            mix = E.GasMixture(air, n)
            alloc_sum(mix, n) = @allocated TASOPT.engine.gassumd(mix, n, 900.0)
            alloc_new(a, n) = @allocated TASOPT.engine.GasMixture(a, n)
            alloc_sum(mix, n); alloc_new(air, n)
            @test alloc_sum(mix, n) == 0
            @test alloc_new(air, n) == 0
            @test_throws ArgumentError E.gassum(mix, n - 1, 900.0)
            @test_throws ArgumentError E.GasMixture(ones(6), 6)
            # derivatives with respect to temperature and composition
            @test TASOPT.ForwardDiff.derivative(x -> E.gassum(mix, n, x)[3], 900.0) ≈ E.gassum(air, n, 900.0)[5] rtol = 1e-13
            dlam(f) = E.gassum(E.GasMixture(products(f), n), n, 1400.0)[3]
            @test TASOPT.ForwardDiff.derivative(dlam, 0.02) ≈
                  TASOPT.ForwardDiff.derivative(f -> E.gassum(products(f), n, 1400.0)[3], 0.02) rtol = 1e-12
        end
    end

    @testset "Temperature derivatives are consistent" begin
        alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        n = length(alpha)
        δ = 1e-2
        for T in (220.0, 450.0, 999.0, 1001.0, 1800.0, 2400.0)
            s, s_t, h, h_t, cp, cp_t, r = E.gassumd(alpha, n, T)
            @test h_t ≈ cp rtol = 1e-14
            @test s_t ≈ cp / T rtol = 1e-14
            # analytic dcp/dT against a central difference of cp
            cp_p = E.gassum(alpha, n, T + δ)[5]
            cp_m = E.gassum(alpha, n, T - δ)[5]
            @test cp_t ≈ (cp_p - cp_m) / (2δ) rtol = 1e-6
            # automatic differentiation through gassum agrees with the
            # analytic derivatives (the legacy table code rejected Dual numbers)
            @test TASOPT.ForwardDiff.derivative(x -> E.gassum(alpha, n, x)[3], T) ≈ cp rtol = 1e-12
            @test TASOPT.ForwardDiff.derivative(x -> E.gassum(alpha, n, x)[5], T) ≈ cp_t rtol = 1e-10
            # dh/dT and ds/dT against central differences of h and s
            @test h_t ≈ (E.gassum(alpha, n, T + δ)[3] - E.gassum(alpha, n, T - δ)[3]) / (2δ) rtol = 1e-8
            @test s_t ≈ (E.gassum(alpha, n, T + δ)[1] - E.gassum(alpha, n, T - δ)[1]) / (2δ) rtol = 1e-8
        end
    end
end
