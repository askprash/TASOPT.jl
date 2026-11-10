
@testset "PEM fuel cell" begin

    u_LT = TASOPT.engine.PEMFC_inputs()

    u_LT.j = 1e4
    u_LT.T = 353.15
    u_LT.p_A = 3e5
    u_LT.p_C = 3e5
    u_LT.x_H2O_A = 0.1
    u_LT.x_H2O_C = 0.1
    u_LT.λ_H2 = 3.0
    u_LT.λ_O2 = 3.0
    u_LT.t_M = 100e-6
    u_LT.t_A = 250e-6
    u_LT.t_C = 250e-6
    u_LT.type = "LT-PEMFC"

    V_LT, α_LT = TASOPT.engine.LT_PEMFC_voltage(u_LT)

    @test V_LT ≈ 0.7103339015901087
    @test α_LT ≈ 0.13026807266756726

    u_HT = TASOPT.engine.PEMFC_inputs()

    u_HT.j = 1e4
    u_HT.T = 453.15
    u_HT.p_A = 3e5
    u_HT.p_C = 3e5
    u_HT.x_H2O_A = 0.1
    u_HT.x_H2O_C = 0.1
    u_HT.λ_H2 = 3.0
    u_HT.λ_O2 = 3.0
    u_HT.t_M = 100e-6
    u_HT.t_A = 250e-6
    u_HT.t_C = 250e-6
    u_HT.type = "HT-PEMFC"

    V_HT= TASOPT.engine.HT_PEMFC_voltage(u_HT)

    @test V_HT ≈ 0.7532410548758377

    Pdes = 1e6
    Vdes = 200

    n_cells, A_cell, Q = TASOPT.engine.PEMsize(Pdes, Vdes, u_HT)

    @test n_cells ≈ 265.5192500533146
    @test A_cell ≈ 0.5
    @test Q ≈ 664805.6978342824

    V_stack, Q_od = TASOPT.engine.PEMoper(5e5, n_cells, A_cell, u_HT)

    @test V_stack ≈ 237.6830551160271
    @test Q_od ≈ 200430.95710865577

    W_FC = TASOPT.engine.PEMstackweight(9.81, u_HT, n_cells, A_cell, 4.0)

    @test W_FC ≈ 7032.808376162144
end

# LT-PEMFC voltage and water-flux ratio, against reference solutions of the membrane ODEs (OrdinaryDiffEq Vern9 at
# reltol = abstol = 1e-13). Cases: (j, T, p, t_M) or (j, T, p, t_M, x_H2O_A, x_H2O_C, λ); the fifth case is
# unphysical, so its voltage is clamped to 0.
@testset "LT-PEMFC membrane ODEs" begin
    function lt_inputs(j, T, p, t_M, x_H2O_A = 0.1, x_H2O_C = 0.1, λ = 3.0)
        u = TASOPT.engine.PEMFC_inputs()
        u.j = j
        u.T = T
        u.p_A = p
        u.p_C = p
        u.x_H2O_A = x_H2O_A
        u.x_H2O_C = x_H2O_C
        u.λ_H2 = λ
        u.λ_O2 = λ
        u.t_M = t_M
        u.t_A = 250e-6
        u.t_C = 250e-6
        u.type = "LT-PEMFC"
        return u
    end
    refs = [
        ((10000.0, 353.15, 300000.0, 0.0001), 0.71033381436589593, 0.13026808012786686),
        ((1000.0, 353.15, 300000.0, 0.0001), 0.95229858971855319, -0.025798977880345606),
        ((5000.0, 343.15, 200000.0, 5.0e-5), 0.85896175593059976, 0.018938558907262856),
        ((20000.0, 353.15, 300000.0, 0.0001), 0.21273352343891316, 0.26205101899927713),
        ((10000.0, 333.15, 150000.0, 0.000175), 0.0, 0.29021159368525046),
        ((3000.0, 363.15, 250000.0, 2.5e-5), 0.89766502445096752, -0.15593556294911934),
        ((15000.0, 353.15, 300000.0, 0.0001, 0.3, 0.05, 2.0), 0.56913245697534554, 0.73405499843619193),
    ]
    for (case, V_ref, α_ref) in refs
        V, α_star = TASOPT.engine.LT_PEMFC_voltage(lt_inputs(case...))
        @test V ≈ V_ref rtol = 1e-5
        @test α_star ≈ α_ref rtol = 1e-5
    end
end
