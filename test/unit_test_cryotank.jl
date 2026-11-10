using StaticArrays
using Random

@testset "Cryogenic tank" begin
    p_atm = TASOPT.p_atm
    Δp = 0.001 * p_atm
    @testset "Thermodynamic properties" begin
        #Hydrogen
        species = "H2"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [0.71653, 1.3322, 2.5133, 6.1715]
        us = [364.55, 372.65, 379.15, 378.1] * 1e3
        hs = [435.25, 448.71, 459.78, 460.2] * 1e3
        
        @testset "Hydrogen gas" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.gas_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.5
                @test props.ρ ≈ ρs[i] atol = 1e-1
                @test props.u ≈ us[i] rtol = 1e-2
                @test props.h ≈ hs[i] rtol = 1e-2

                #Check derivative by finite differences
                props2 = TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [73.142, 70.848, 67.639, 60.711]
        us = [-20.876, -1.4302, 25.31, 80.495] * 1e3
        hs = [-20.183, 0, 28.306, 88.84] * 1e3

        @testset "Hydrogen liquid" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.liquid_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.5
                @test props.ρ ≈ ρs[i] rtol = 1e-2
                @test props.u ≈ us[i] atol = 2e3
                @test props.h ≈ hs[i] atol = 2e3

                props2 = TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Methane
        species = "CH4"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [0.96233, 1.8164, 3.4382, 8.1024]
        us = [444.65, 455.05, 466.28, 481.52] * 1e3
        hs = [497.3, 510.83, 525.22, 544.05] * 1e3

        @testset "Methane gas" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.gas_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.6
                @test props.ρ ≈ ρs[i] atol = 1e-1
                @test props.u ≈ us[i] rtol = 1e-2
                @test props.h ≈ hs[i] rtol = 1e-2

                props2 = TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [433.5, 422.36, 408.66, 384.63]
        us = [-27.112, 0, 31.803, 85.119] * 1e3
        hs = [-26.995, 0, 32.299, 86.436] * 1e3

        @testset "Methane liquid" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.liquid_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.6
                @test props.ρ ≈ ρs[i] rtol = 1e-2
                @test props.u ≈ us[i] atol = 2e3
                @test props.h ≈ hs[i] atol = 2e3

                props2 = TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end
    end

    @testset "Saturated mixtures" begin
        species = "H2"
        
        p = p_atm
        β = 0.5

        @testset "SaturatedMixture assembly" begin
            # gas_properties / liquid_properties → SaturatedPhaseProps stored directly
            # in mix.gas and mix.liquid; verify the constructor assembles them correctly.
            gas = TASOPT.CryoTank.gas_properties(species, p)
            liq = TASOPT.CryoTank.liquid_properties(species, p)
            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            @test mix.gas === mix.gas       # sanity: field exists
            @test mix.p ≈ p
            @test mix.T ≈ gas.Tsat         # temperature comes from saturation curve

            # Verify mix.gas and mix.liquid hold the phase properties
            @test mix.gas.Tsat ≈ gas.Tsat
            @test mix.gas.ρ    ≈ gas.ρ
            @test mix.gas.h    ≈ gas.h
            @test mix.gas.u    ≈ gas.u
            @test mix.gas.ρ_p  ≈ gas.ρ_p
            @test mix.gas.u_p  ≈ gas.u_p
            @test mix.liquid.ρ    ≈ liq.ρ
            @test mix.liquid.h    ≈ liq.h
            @test mix.liquid.u    ≈ liq.u
            @test mix.liquid.ρ_p  ≈ liq.ρ_p
            @test mix.liquid.u_p  ≈ liq.u_p

            # Verify the mixture-level scalars against the mixing rule
            x = 1 / (1 + (liq.ρ/gas.ρ) * (β / (1 - β)))
            ρ = 1 / (x / gas.ρ + (1-x) / liq.ρ)
            u = x * gas.u + (1-x) * liq.u
            h = x * gas.h + (1-x) * liq.h
            x_p = (-x/gas.ρ^2 * gas.ρ_p - (1-x)/liq.ρ^2 * liq.ρ_p) / (1/liq.ρ - 1/gas.ρ)
            u_p = x * gas.u_p + (1 - x) * liq.u_p + x_p * (gas.u - liq.u)
            ϕ      = 1 / (u_p * ρ)
            hvap   = gas.h - liq.h
            ρ_star = gas.ρ / (liq.ρ - gas.ρ)

            @test mix.ρ     ≈ ρ
            @test mix.x     ≈ x
            @test mix.h     ≈ h
            @test mix.u     ≈ u
            @test mix.β     ≈ β
            @test mix.u_p   ≈ u_p
            @test mix.ϕ     ≈ ϕ
            @test mix.hvap  ≈ hvap
            @test mix.ρ_star ≈ ρ_star
        end
        @testset "Mixture functions" begin     

            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            #Test update to mixture by comparing against new object
            p2 = 2 * p_atm 
            β2 = 0.3
            TASOPT.CryoTank.update_pβ!(mix, p2, β2)

            mix2 = TASOPT.CryoTank.SaturatedMixture(species, p2, β2)

            @test mix2.p ≈ mix.p
            @test mix2.T ≈ mix.T
            @test mix2.ρ ≈ mix.ρ
            @test mix2.x ≈ mix.x
            @test mix2.h ≈ mix.h
            @test mix2.u ≈ mix.u
            @test mix2.β ≈ mix.β
            @test mix2.u_p ≈ mix.u_p
            @test mix2.ϕ ≈ mix.ϕ
            @test mix2.hvap ≈ mix.hvap
            @test mix2.ρ_star ≈ mix.ρ_star

            #Check that density is the same during pressure change with convert_β_same_ρ
            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)
            ρorig = mix.ρ

            β3 = TASOPT.CryoTank.convert_β_same_ρ(species, p2, p, β)
            TASOPT.CryoTank.update_pβ!(mix, p2, β3)
            @test mix.ρ ≈ ρorig

        end
    end

    species = "h2"
    p = 2.0 * p_atm
    β = 0.5 
    mix = TASOPT.CryoTank.SaturatedMixture("h2", p, β)

    Q = 2e4
    mdot = 1.0
    W = 0.0
    mdot_vent = 0.0

    xout = 0.0 #Taking liquid out
    xvent = 1.0 #Venting gas
    α = 2.0
    V = 100.0
    pmax = 5*p_atm

    @testset "Model derivatives" begin
        #Check function for dp/dt
        dp_dt = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent, xvent, V, α)
        dp_dt_validate = α * mix.ϕ / V * (Q + W - mdot * mix.hvap * (xout + mix.ρ_star) - mdot_vent * mix.hvap * (xvent + mix.ρ_star))
        dp_dt_regression = 6.995703982056127

        @test dp_dt ≈ dp_dt_validate
        @test dp_dt ≈ dp_dt_regression

        #Check function for dβ/dt
        dβ_dt = TASOPT.CryoTank.dβdt(mix, dp_dt_regression, mdot + mdot_vent, V)
        dβ_dt_validate = (-(mdot + mdot_vent) / V - dp_dt_validate * (β * mix.liquid.ρ_p + (1-β)*mix.gas.ρ_p)) / (mix.liquid.ρ - mix.gas.ρ)
        dβ_dt_regression = -0.00015294953663597238

        @test dβ_dt ≈ dβ_dt_validate
        @test dβ_dt ≈ dβ_dt_regression

        #Check function for venting mass flow rate
        mdot_vent_out = TASOPT.CryoTank.venting_mass_flow(mix, Q, W, mdot, xout, xvent) #This venting should give zero dp/dt
        vent_validate = (Q - mdot * mix.hvap * (xout + mix.ρ_star)) / (mix.hvap * (xvent + mix.ρ_star))
        vent_regression = 0.007274844258979718

        @test mdot_vent_out ≈ vent_validate
        @test mdot_vent_out ≈ vent_regression

        dpdt_venting = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent_out, xvent, V, α) #Check pressure derivative with venting; it should be 0
        dpdt_venting_sanity = 0.0

        @test dpdt_venting ≈ dpdt_venting_sanity

        #Check function for boiloff rate
        mdot_liq = (1 - xout) * mdot + (1 - xvent) * mdot_vent
        mdotboil = TASOPT.CryoTank.mdot_boiloff(mix, dβ_dt_regression, dp_dt_regression, mdot_liq, V)
        mdotboil_validate = - (dβ_dt_regression * V * mix.liquid.ρ + β * V * mix.liquid.ρ_p * dp_dt_regression + mdot)
        mdotboil_regression = 0.04266684331352266

        @test mdotboil ≈ mdotboil_validate
        @test mdotboil ≈ mdotboil_regression

        #Check function that computes all derivatives for inputs
        Q_calc(t) = Q
        W_calc(t) = W
        mdot_calc(t) = mdot

        u = TASOPT.CryoTank.tank_inputs(Q_calc, W_calc, mdot_calc)
        params = TASOPT.CryoTank.tank_params(mix, V, pmax, xout, xvent, α)

        y = @SVector [p, β, V * mix.ρ, 0.0, 0.0, 0.0]
        t = 0.0

        dydt = TASOPT.CryoTank.TankDerivatives(t, y, u, params)
        @test dydt[1] ≈ dp_dt_regression
        @test dydt[2] ≈ dβ_dt_regression
        @test dydt[3] ≈ -mdot
        @test dydt[4] ≈ mdot
        @test dydt[5] ≈ 0.0
        @test dydt[6] ≈ mdotboil_regression

        #The model is not valid once the tank has no liquid left: it throws, even at the vent pressure
        params_vent = TASOPT.CryoTank.tank_params(mix, V, p, xout, xvent, α) #pmax = p: the tank vents
        for β_empty in (0.0, -1e-3)
            y_empty = @SVector [p, β_empty, 10.0, 1e3, 50.0, 20.0]
            @test_throws ErrorException TASOPT.CryoTank.TankDerivatives(t, y_empty, u, params_vent)
        end
        #Just before it empties, the tank still drains and vents
        dydt_full = TASOPT.CryoTank.TankDerivatives(t, (@SVector [p, 1e-3, 10.0, 1e3, 50.0, 20.0]), u, params_vent)
        @test dydt_full[4] ≈ mdot
        @test dydt_full[5] > 0.0
    end
end

# Tank state over a mission, against reference solutions (OrdinaryDiffEq Vern9 at reltol = abstol = 1e-13). The inputs
# are synthetic but built like analyze_TASOPT_tank builds them: mission points with leading zeros, a negative time
# outside the span, duplicated segment boundaries (jumps in the fuel flow) and a hold at origin. The two H2 missions
# draw more fuel than the tank holds: the first runs dry in cruise without venting, the second vents during a two-hour
# hold and then runs dry. Integrating them to the end must throw; up to the last mission point before dry-out they
# are checked against the references. The methane case neither vents nor runs dry.
@testset "Tank state integration" begin
    CT = TASOPT.CryoTank
    function tank_case(species, hold, Qscale)
        tm = [0.0, 0.0, -30.0, 0.0, 0.0, 0.0, 0.0, 180.0, 360.0, 630.0, 1200.0, 1200.0, 22800.0, 22800.0,
              23000.0, 23250.0, 23500.0, 24000.0, 24000.0, 24000.0]
        md = [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.18, 0.15, 0.12, 0.1, 0.07, 0.06, 0.04, 0.035, 0.03, 0.03,
              0.025, 0.0, 0.0]
        Qp = Qscale .* [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.95, 0.9, 0.85, 0.8, 0.8, 0.8, 0.8, 0.85, 0.9, 0.95,
                        1.0, 1.0, 1.0]
        times = copy(tm)
        times[2:end] .+= hold
        times[1] = 0.0
        para = zeros(iatotal, length(times))
        para[iatime, :] .= times
        pare = zeros(ietotal, length(times))
        pare[ieff, :] .= 1.0
        pare[iemcore, :] .= md
        parg = zeros(igtotal)
        parg[igneng] = 2.0
        Q_calc(t::Float64) = CT.find_Q_time_interp(t, para, Qp)
        W_calc(t::Float64) = 0.0
        mdot_calc(t::Float64) = CT.find_mdot_time(t, 1, parg, para, pare)
        u = CT.tank_inputs(Q_calc, W_calc, mdot_calc)
        V = 30.0
        pmax = 2.0e5
        p0 = 1.2e5
        β0 = 0.9
        mix = CT.SaturatedMixture(species, p0, β0)
        params = CT.tank_params(mix, V, pmax, 0.0, 1.0, 1.0)
        y0 = SVector(p0, β0, V * mix.ρ, 0.0, 0.0, 0.0)
        return u, params, y0, (0.0, times[end]), times
    end

    # States [p, β, M, Mburn, Mvent, Mboil] at every mission point inside the span
    refs = Dict(
        ("H2", 0.0, 2000.0) => [
            (180.0, [119435.28373403501, 0.86653104672332482, 1829.6231223806451, 68.336795383414696, 0.0, 1.5238409689898378]),
            (360.0, [118987.98906603278, 0.83755876779302818, 1770.3871209449767, 127.57279681908314, 0.0, 2.8344117110433884]),
            (630.0, [118544.67625327989, 0.80211784923788587, 1697.7881150378391, 200.17180272622085, 0.0, 4.4323141153086691]),
            (1200.0, [117993.03077622077, 0.74119861449369673, 1572.7343342292061, 325.22558353485414, 0.0, 7.1721166542361514]),
        ],
        ("H2", 7200.0, 20000.0) => [
            (7170.0, [200000.00000000114, 0.81420296753256194, 1666.0065801175861, 0.0, 231.95333764647225, 241.20168900058792]),
            (7200.0, [200000.00000000114, 0.81351886505522686, 1664.669560606713, 0.0, 233.29035715734534, 242.5898224867741]),
            (7380.0, [200000.00000000114, 0.77583900138033768, 1591.0275046729, 68.336795383414, 238.59561770774468, 250.71040338096321]),
            (7560.0, [200000.00000000114, 0.7428494383265275, 1526.5522382918336, 127.5727968190825, 243.83488265314256, 258.41454422929655]),
            (7830.0, [200000.00000000114, 0.70168375676825878, 1446.0974508216805, 200.17180272622016, 251.69066421615807, 269.34609543301269]),
            (8400.0, [200000.00000000114, 0.62933121801520819, 1304.6906306755566, 325.22558353485357, 268.04370355364853, 291.10508782898108]),
        ],
        ("CH4", 1800.0, 5000.0) => [
            (1770.0, [121973.20161001217, 0.9007011916333435, 11325.545503843301, 0.0, 0.0, 0.050598739028490353]),
            (1800.0, [122006.88542219321, 0.90071309116463483, 11325.545503843301, 0.0, 0.0, 0.051454928122483888]),
            (1980.0, [122164.63203318696, 0.89530269799049234, 11257.208708459892, 68.336795383414653, 0.0, 0.40625283521056965]),
            (2160.0, [122318.52870580571, 0.89061816328763999, 11197.972707024224, 127.57279681908298, 0.0, 0.71503466467972121]),
            (2430.0, [122545.17302371706, 0.88488904691183778, 11125.373701117087, 200.17180272622053, 0.0, 1.0957593210157812]),
            (3000.0, [123012.7837332482, 0.87504381648750496, 11000.319920308455, 325.22558353485391, 0.0, 1.7576757050080445]),
            (24600.0, [145665.2598500218, 0.65594560846354488, 8197.8671482277405, 3127.6783556155719, 0.0, 19.668483543090694]),
            (24800.0, [145944.24419024211, 0.65479513765408115, 8182.8893968474267, 3142.6561069958889, 0.0, 19.801685873497547]),
            (25050.0, [146319.51405222909, 0.65356443857184343, 8166.6714988608492, 3158.8740049824655, 0.0, 19.957316624905086]),
            (25300.0, [146720.24445034305, 0.65243693207723341, 8151.6714988608492, 3173.8740049824655, 0.0, 20.109797470783874]),
            (25800.0, [147575.09034008649, 0.65039951718036892, 8124.2474241221144, 3201.2980797212008, 0.0, 20.408389323418731]),
        ],
    )

    runs_dry = Dict(("H2", 0.0, 2000.0) => true, ("H2", 7200.0, 20000.0) => true, ("CH4", 1800.0, 5000.0) => false)

    for ((species, hold, Qscale), yref) in refs
        u, params, y0, tspan, tstops = tank_case(species, hold, Qscale)
        if runs_dry[(species, hold, Qscale)]
            @test_throws ErrorException CT.integrate_tank_state(u, params, y0, tspan, tstops)
            tspan = (tspan[1], yref[end][1]) #up to the last mission point before the tank runs dry
        end
        ts, ys = CT.integrate_tank_state(u, params, y0, tspan, tstops)
        @test ts[1] == tspan[1]
        @test ts[end] == tspan[2]
        @test issorted(ts)
        # Errors are measured against the scale of each state (initial pressure, 1, initial mass)
        scale = SVector(y0[1], 1.0, y0[3], y0[3], y0[3], y0[3])
        for (t, y) in yref
            i = findfirst(==(t), ts) # the integration must step onto every mission point
            @test i !== nothing
            i === nothing && continue
            @test all(abs.(ys[i] .- y) .<= 1e-5 .* scale)
        end
    end
end

# ── Property-based tests ────────────────────────────────────────────────────
# Inspired by https://fsharpforfunandprofit.com/posts/property-based-testing-2/
# Strategy: test *relationships* and *invariants* that must hold for all valid
# inputs.
@testset "SaturatedPhaseProps property tests" begin
    rng = MersenneTwister(42)
    p_atm_val = TASOPT.p_atm
    p_lo = 0.1 * p_atm_val   # lower bound of polynomial fit domain
    p_hi = 10.0 * p_atm_val  # upper bound

    for species in ("H2", "CH4")
        @testset "$(species): saturation temperature" begin
            # gas_properties and liquid_properties are both fits to the same
            # saturation curve - Tsat must be identical at any given pressure.
            for _ in 1:10
                p = p_lo + rand(rng) * (p_hi - p_lo)
                gp = TASOPT.CryoTank.gas_properties(species, p)
                lp = TASOPT.CryoTank.liquid_properties(species, p)
                @test gp.Tsat ≈ lp.Tsat rtol = 1e-12
            end
        end

        @testset "$(species): Tsat monotone with pressure (Clausius-Clapeyron)" begin
            for _ in 1:10
                p1 = p_lo + rand(rng) * (p_hi - p_lo - 0.1 * p_atm_val)
                p2 = p1 + 0.05 * p_atm_val + rand(rng) * 0.5 * p_atm_val
                gp1 = TASOPT.CryoTank.gas_properties(species, p1)
                gp2 = TASOPT.CryoTank.gas_properties(species, p2)
                @test gp2.Tsat > gp1.Tsat
            end
        end

        @testset "$(species): thermodynamic identity h = u + p/ρ" begin
            # Follows from the definition h ≡ u + pv, v = 1/ρ.
            # h and u are independent NIST polynomial fits, so the identity holds
            # to ~1% inside the fit range (0.5–9 atm) but degrades at the edges.
            # We stay away from the edges.
            p_lo_inner = 0.5 * p_atm_val
            p_hi_inner = 9.0 * p_atm_val
            for _ in 1:10
                p = p_lo_inner + rand(rng) * (p_hi_inner - p_lo_inner)
                for props in (TASOPT.CryoTank.gas_properties(species, p),
                              TASOPT.CryoTank.liquid_properties(species, p))
                    @test props.h ≈ props.u + p / props.ρ rtol = 1e-2
                end
            end
        end

        @testset "$(species): derivative consistency" begin
            # ρ_p and u_p are analytic derivatives of the same polynomials.
            # Try one-sided finite difference to within 0.05 %.
            δ = 10.0  # 10 Pa — small relative to ~1e5 Pa but large enough to avoid cancellation
            for _ in 1:10
                p = p_lo + rand(rng) * (p_hi - p_lo - δ)
                for (f, label) in ((TASOPT.CryoTank.gas_properties, "gas"),
                                   (TASOPT.CryoTank.liquid_properties, "liquid"))
                    props  = f(species, p)
                    props2 = f(species, p + δ)
                    @test props.ρ_p ≈ (props2.ρ - props.ρ) / δ rtol = 5e-4
                    @test props.u_p ≈ (props2.u - props.u) / δ rtol = 5e-4
                end
            end
        end
    end

    
end