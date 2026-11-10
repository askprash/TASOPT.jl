using StaticArrays

@testset "ode_solve" begin
    ode_solve = TASOPT.ode_solve

    @testset "Scalar decay" begin
        # dy/dt = -y, y(0) = 1: y = exp(-t)
        for (rtol, err) in ((1e-6, 1e-5), (1e-10, 1e-9))
            ts, ys, ok = ode_solve((y, t) -> -y, 1.0, (0.0, 5.0); reltol = rtol, abstol = rtol)
            @test ok
            @test ts[1] == 0.0 && ts[end] == 5.0
            @test maximum(abs.(ys .- exp.(-ts))) < err
        end
    end

    @testset "Harmonic oscillator, SVector state" begin
        # y1' = y2, y2' = -y1, y(0) = (1, 0): y = (cos t, -sin t)
        f(y, t) = SVector(y[2], -y[1])
        ts, ys, ok = @inferred ode_solve(f, SVector(1.0, 0.0), (0.0, 10.0); reltol = 1e-10, abstol = 1e-10)
        @test ok
        @test eltype(ys) == SVector{2, Float64}
        @test maximum(maximum(abs.(y .- SVector(cos(t), -sin(t)))) for (t, y) in zip(ts, ys)) < 1e-8
        # Fewer steps at a looser tolerance
        ts6, _, _ = ode_solve(f, SVector(1.0, 0.0), (0.0, 10.0))
        @test length(ts6) < length(ts)
    end

    @testset "Time-dependent derivative" begin
        # dy/dt = cos(t), y(0) = 0: y = sin(t)
        ts, ys, ok = ode_solve((y, t) -> cos(t), 0.0, (0.0, 20.0); reltol = 1e-10, abstol = 1e-10)
        @test ok
        @test maximum(abs.(ys .- sin.(ts))) < 1e-8
    end

    @testset "tstops" begin
        # Piecewise-constant derivative with jumps at the stops: every step is exact, whichever value f takes
        # exactly at a jump
        f(y, t) = t < 1.0 ? 1.0 : (t < 2.5 ? -2.0 : 0.5)
        g(y, t) = t <= 1.0 ? 1.0 : (t <= 2.5 ? -2.0 : 0.5)
        tsg, ysg, okg = ode_solve(g, 0.0, (0.0, 3.0); tstops = [1.0, 2.5])
        @test okg
        @test ysg[end] ≈ -1.75 atol = 1e-12
        # A jump at the initial time
        tsj, ysj, okj = ode_solve((y, t) -> t <= 0.0 ? 100.0 : 1.0, 0.0, (0.0, 1.0))
        @test okj
        @test ysj[end] ≈ 1.0 atol = 1e-12
        tstops = [2.5, -1.0, 1.0, 1.0, 0.0, 3.0, 7.0] # unsorted, repeated, at and outside the ends
        ts, ys, ok = ode_solve(f, 0.0, (0.0, 3.0); tstops = tstops)
        @test ok
        @test issorted(ts) && allunique(ts)
        @test ts[1] == 0.0 && ts[end] == 3.0
        @test 1.0 in ts && 2.5 in ts
        @test ys[findfirst(==(1.0), ts)] ≈ 1.0 atol = 1e-12
        @test ys[findfirst(==(2.5), ts)] ≈ -2.0 atol = 1e-12
        @test ys[end] ≈ -1.75 atol = 1e-12
    end

    @testset "Edge cases" begin
        ts, ys, ok = ode_solve((y, t) -> -y, 2.0, (1.0, 1.0))
        @test ok && ts == [1.0] && ys == [2.0]
        @test_throws ArgumentError ode_solve((y, t) -> -y, 1.0, (1.0, 0.0))
        # Solution blows up at t = 1 (y = 1/(1 - t)): the integration stops with a warning
        ts, ys, ok = @test_logs (:warn,) match_mode = :any ode_solve((y, t) -> y^2, 1.0, (0.0, 2.0))
        @test !ok
        @test ts[end] < 2.0
    end
end
