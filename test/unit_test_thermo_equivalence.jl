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
            # dh/dT and ds/dT against central differences of h and s
            @test h_t ≈ (E.gassum(alpha, n, T + δ)[3] - E.gassum(alpha, n, T - δ)[3]) / (2δ) rtol = 1e-8
            @test s_t ≈ (E.gassum(alpha, n, T + δ)[1] - E.gassum(alpha, n, T - δ)[1]) / (2δ) rtol = 1e-8
        end
    end
end
