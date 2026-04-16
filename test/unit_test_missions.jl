@testset "mission" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    #a. verify that fly_mission!() produces same results as size_aircraft!() call
    pfei_sizecall = ac.parm[imPFEI, 1]

    fly_mission!(ac, 1; printTO=false)
    pfei_missioncall = ac.parm[imPFEI, 1]
    @test pfei_sizecall ≈ pfei_missioncall

end

@testset "aircraft missions field" begin
    ac = load_default_model()
    # missions is populated with one entry per N_missions (default input has 2)
    nmisx = size(ac.pare, 3)
    @test length(ac.missions) == nmisx
    # Each mission has iptotal (17) points
    @test length(ac.missions[1].points) == 17
    # design_mission_state shortcut returns the design Mission
    @test ac.design_mission_state isa Mission{Float64}
    @test ac.design_mission_state === ac.missions[1]
end

@testset "MissionPoint and Mission constructors" begin
    # Mission(npoints) constructs without error and has correct length
    m = Mission(17)
    @test length(m.points) == 17

    # Each point holds a Float64 EngineState
    @test m.points[1].engine isa TASOPT.engine.EngineState{Float64}

    # Parametric constructor
    m32 = Mission{Float32}(3)
    @test length(m32.points) == 3
    @test m32.points[1].engine isa TASOPT.engine.EngineState{Float32}

    # Empty constructor
    m0 = Mission()
    @test length(m0.points) == 0

    # MissionPoint default constructor
    pt = MissionPoint()
    @test pt.engine isa TASOPT.engine.EngineState{Float64}
end
