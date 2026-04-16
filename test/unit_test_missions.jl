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

@testset "fly_mission! mirrors pare into ac.missions engine state" begin
    # After fly_mission! each mission point that goes through _mission_iteration!'s
    # enginecalc! must have its ac.missions[im].points[ip].engine synchronised with
    # the corresponding pare[:, ip, im] column.
    ac = load_default_model()
    size_aircraft!(ac; printiter=false)
    fly_mission!(ac, 1; printTO=false)

    im = 1
    tol = 1e-12

    # ---- ipcruise1: three canonical checks from the issue spec ----
    eng_cr = ac.missions[im].points[ipcruise1].engine
    @test eng_cr.TSFC ≈ ac.pare[ieTSFC, ipcruise1, im]  rtol=tol
    @test eng_cr.Fe   ≈ ac.pare[ieFe,   ipcruise1, im]  rtol=tol
    @test eng_cr.st4.Tt ≈ ac.pare[ieTt4, ipcruise1, im] rtol=tol

    # ---- ipcruisen: end-of-cruise point ----
    eng_cr2 = ac.missions[im].points[ipcruisen].engine
    @test eng_cr2.TSFC ≈ ac.pare[ieTSFC, ipcruisen, im] rtol=tol
    @test eng_cr2.Fe   ≈ ac.pare[ieFe,   ipcruisen, im] rtol=tol

    # ---- climb segment: all points ipclimb1:ipclimbn ----
    for ip in ipclimb1:ipclimbn
        eng_ip = ac.missions[im].points[ip].engine
        @test eng_ip.TSFC ≈ ac.pare[ieTSFC, ip, im] rtol=tol
        @test eng_ip.Fe   ≈ ac.pare[ieFe,   ip, im] rtol=tol
    end

    # ---- descent segment: all points ipdescent1:ipdescentn ----
    for ip in ipdescent1:ipdescentn
        eng_ip = ac.missions[im].points[ip].engine
        @test eng_ip.TSFC ≈ ac.pare[ieTSFC, ip, im] rtol=tol
        @test eng_ip.Fe   ≈ ac.pare[ieFe,   ip, im] rtol=tol
    end

    # ---- freestream p0 at ipcruisen (read via typed state in _mission_iteration!) ----
    eng_cn = ac.missions[im].points[ipcruisen].engine
    @test eng_cn.p0 ≈ ac.pare[iep0, ipcruisen, im] rtol=tol
    # mfuel at covered call sites
    @test ac.missions[im].points[ipcruisen].engine.mfuel ≈ ac.pare[iemfuel, ipcruisen, im] rtol=tol
    for ip in ipclimb1:ipclimbn
        @test ac.missions[im].points[ip].engine.mfuel ≈ ac.pare[iemfuel, ip, im] rtol=tol
    end
    for ip in ipdescent1:ipdescentn
        @test ac.missions[im].points[ip].engine.mfuel ≈ ac.pare[iemfuel, ip, im] rtol=tol
    end
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
