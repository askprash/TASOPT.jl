

import ForwardDiff

@testset "wing aerodynamics" verbose=true begin
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))
    fuselage = ac_test.fuselage
    wing = ac_test.wing
    airf = wing.airsection
    @testset "Airfoil data" begin
        #Test that it is reading the database correctly first
        @test airf.Re == 2.0e7
        @test airf.cl[1] == 0.4
        @test airf.cl[end] == 0.9
        @test airf.τ[1] == 0.09
        @test airf.τ[end] == 0.145
        @test airf.Ma[1] == 0.0
        @test airf.Ma[end] == 0.81
        @test all(airf.A[1, 1, 1, :] .== [0.537000E-02, 0.105000E-02, -0.114800])
        @test all(airf.A[end, end, end, :] .== [0.330800E-02, 0.129780, -0.157000E-01])

        # Airfun:
        clp = 0.68657968661106417
        toc = 0.12677500000000000
        Mperp = 0.73161835368752204
        cdf1 = 4.9548669970280916E-003
        cdp1 = 3.8227447198802837E-003
        cdwbar = 0.0000000000000000
        cm1 = -0.11614728079885453

        cdf, cdp, cdw, cm = TASOPT.aerodynamics.airfun(clp, toc, Mperp, airf)

        @test cdf1 ≈ cdf rtol = 1e-6
        @test cdp1 ≈ cdp rtol = 1e-6
        @test cdwbar ≈ cdw
        @test cm1 ≈ cm rtol = 1e-6

        clp = 0.67016181769975336
        toc = 0.14401500000000000
        Mperp = 0.72079998136376022
        cdf1 = 4.9388801313136523E-003
        cdp1 = 4.1045733512642419E-003
        cdwbar = 0.0000000000000000
        cm1 = -9.5253127678246632E-002
        cdf, cdp, cdw, cm = TASOPT.aerodynamics.airfun(clp, toc, Mperp, airf)

        @test cdf1 ≈ cdf
        @test cdp1 ≈ cdp
        @test cdwbar ≈ cdw
        @test cm1 ≈ cm

        @test all(TASOPT.aerodynamics.airfun(0.4, 0.09, 0.3, airf) .== (0.00533, 0.0011, 0.0, -0.1171))
        @test all(TASOPT.aerodynamics.airfun(0.4, 0.09, 0.8, airf) .== (0.00487, 0.00722, 0.0, -0.1544))
        # test_limits
        cdf, cdp, cdw, cm = TASOPT.aerodynamics.airfun(0.0, 0.0, 0.8, airf)

        @test cdp > 0.0

    end # end airfun

    #start wing_loading
    wing.layout.root_chord = 5.3938688126436549
    wing.layout.span = 35.486921629195265
    wing.layout.ηs = 10.113772664320649 / wing.layout.span
    wing.layout.root_span = 3.6067999999999998
    wing.layout.S = 105.24594887917868
    wing.inboard.λ = 0.69999999999999996
    wing.outboard.λ = 0.25000000000000000
    wing.layout.AR = 10.1
    wing.fuse_lift_carryover = -0.29999999999999998890 
    wing.tip_lift_loss = -0.05000000000000000278 
    
    rclt = 0.9
    rcls = 1.238
    N = 3.0000000000000000
    W = 778345.75427325454  
    Lhtail = -132476.65894384126  
    fort_po = 110091.58394892939 

    po = TASOPT.aerodynamics.wing_loading(wing, rclt, rcls, N, W, Lhtail)

    @test po ≈ fort_po
    #end wing_loading
    # wing_profiledrag_direct
    #Start wing_profiledrag_direct

    wing.layout.root_chord = 5.3938688126436549
    wing.layout.span = 32.603436685105834
    wing.layout.ηs = 9.2919794552551611 / wing.layout.span
    wing.layout.root_span = 3.6067999999999998
    wing.layout.S = 105.24594887917868
    wing.inboard.λ = 0.69999999999999996
    wing.outboard.λ = 0.25000000000000000

    γt = 0.225
    γs = 0.8665999999999999
    Mach = 0.80000000000000004
    CL = 0.56999999999999995
    CLhtail = -2.1202672443102342E-002
    Reco = 34058000.862060823
    aRexp = -0.14999999999999999
    rkSunsw = 0.50000000000000000
    fexcdw = 1.0200000000000000
    fduo = 1.7999999999999999E-002
    fdus = 1.4000000000000000E-002
    fdut = 4.4999999999999997E-003

    fort_clpo = 0.67137920779128768
    fort_clps = 0.83773792743334130
    fort_clpt = 0.62059185385188131
    fort_cdfw = 5.1945675060967979E-003
    fort_cdpw = 3.5732856303350158E-003
    fort_CDwing = 8.7678531364318128E-003
    fort_CDover = 0.0000000000000000

    wing_drag_components = TASOPT.aerodynamics.wing_profiledrag_direct(
      wing, γt, γs,
      Mach, CL, CLhtail, 
      Reco, aRexp, rkSunsw, fexcdw,
      fduo, fdus, fdut)

#       wing_profiledrag_direct INPUT
# (0.225, 0.8665999999999999, 0.5917830310706261, 0.285, -0.0016806695060863123, 6.5275125903133705e7, -0.15, 0.5, 1.02, 0.018, 0.014, 0.0045)
# wing_profiledrag_direct OUTPUT
# (0.005092435287160669, 0.002484528306027699, 0.007576963593188367, 0.0)
    @test fort_clpo ≈ wing_drag_components.clpo
    @test fort_clps ≈ wing_drag_components.clps
    @test fort_clpt ≈ wing_drag_components.clpt
    @test fort_cdfw ≈ wing_drag_components.CDfwing
    @test fort_cdpw ≈ wing_drag_components.CDpwing
    @test fort_CDwing ≈ wing_drag_components.CDwing
    @test fort_CDover ≈ wing_drag_components.CDover
    #end wing_profiledrag_direct

    #start wing_profiledrag_scaled
    S = 124.68530759570760
    b = 35.486921629195265
    bs = 10.113772664320649
    bo = 3.6067999999999998
    λt = 0.25000000000000000
    λs = 0.69999999999999996
    sweep = 26.000000000000000
    co = 5.8841656099573720
    cdfw = 8.6700000000000006E-003
    cdpw = 3.5700000000000003E-003
    Reco = 36429065.140686862
    Rerefw = 20000000.000000000
    aRexp = -0.14999999999999999
    rkSunsw = 0.50000000000000000
    fCDwcen = 0.0000000000000000
    CDwing = 9.4350192385200850E-003
    CDover = 0.0000000000000000
    
    CDw, CDo = TASOPT.aerodynamics.wing_profiledrag_scaled(S,
    b, bs, bo, λt, λs, sweep, co,
    cdfw, cdpw, Reco, Rerefw, aRexp, rkSunsw,
    fCDwcen)

    @test CDwing == CDw
    @test CDover == CDover

    #end wing_profiledrag_scaled


    #start wingsc
    BW = 853967.1303861982
    CL = 0.56999999999999995
    qinf = 10717.328761811295
    # AR = 10.100000000000000
    # etas = 0.28499999999999998
    # bo = 3.6067999999999998
    # lambdat = 0.25000000000000000
    # lambdas = 0.69999999999999996
    fort_S = 139.79117197415886
    fort_b = 37.57513588716619
    fort_bs = 10.708913727842363
    fort_co = 6.239270087906468

    TASOPT.aerodynamics.set_wing_geometry!(BW,CL,qinf,wing)
    # wingsc(BW,CL,qinf,AR,
    # etas,bo,lambdat,lambdas)
    @test fort_S  == wing.layout.S 
    @test fort_b  == wing.layout.span
    @test fort_bs == wing.layout.break_span
    @test fort_co == wing.layout.root_chord
    #end wingsc
    
    #wing_section_cls
    gammat = 0.15
    gammas = 0.77
    CL,CLhtail = 1.2622355275981707,0.019038222769452273
    fduo,fdus,fdut = 0.018,0.014,0.0045
    
    test_clpo = 1.5734696976792315
    test_clps = 1.7444989594158167
    test_clpt = 0.9696283564047791
    clpo, clps, clpt = TASOPT.aerodynamics.wing_section_cls(wing,gammat,gammas,
                              CL,CLhtail,
	                        fduo,fdus,fdut)
    
    @test test_clpo == clpo
    @test test_clps == clps
    @test test_clpt == clpt
    #end wing_section_cls

    #wing_CM
    (b, bs, bo, sweep, Xaxis,
        λt, λs, γt, γs,
        AR, fLo, fLt, cmpo, cmps, cmpt) = 35.723608571355676, 10.717082571406703, 5.7404, 27.567, 0.4, 0.1503, 0.8784, 0.09018, 0.96624, 10.4411, -0.3, -0.05, -0.2, -0.2, -0.02

    CM0, CM1 = TASOPT.aerodynamics.wing_CM(b, bs, bo, sweep, Xaxis,
        λt, λs, γt, γs,
        AR, fLo, fLt, cmpo, cmps, cmpt)

    @test CM0 == -0.05330921996771545
    @test CM1 == -0.3480891589464312
    #end wing_CM

end

@testset "fuse aerodynamics" begin
    #_axisymm_flow
    xnose = 0.0000000000000000
    xend = 37.795200000000001
    xblend1 = 6.0960000000000001
    xblend2 = 29.565600000000000
    Sfuse = 13.507394174276255
    anose = 1.6499999999999999
    btail = 2.0000000000000000
    ifclose = 0
    Mach = 0.80000000000000004
    nc = 30
    nbldim = 60

    nbl = 47
    iblte = 31
    xbl = [0.0000000000000000 0.10352303008851868 0.41295790037283364 0.92491437766070661 1.6337833596531859 2.5317983294431912 3.6091204470999934 4.8539463460583949 6.2526374532728406 7.7898694162777602 9.4487999999999985 11.211253613810754 13.057920447099994 14.968568031602315 16.922262912553194 18.897599999999997 20.872937087446807 22.826631968397688 24.737279552900009 26.583946386189243 28.346399999999996 30.005330583722237 31.542562546727151 32.941253653941608 34.186079552900011 35.263401670556810 36.161416640346815 36.870285622339296 37.382242099627170 37.691676969911484 37.795200000000001 37.898723030088519 38.208157900372832 38.720114377660707 39.428983359653188 40.326998329443192 41.404320447099991 42.649146346058394 44.047837453272848 45.585069416277761 47.244000000000007 49.006453613810763 50.853120447099997 52.763768031602311 54.717462912553195 56.692800000000005 58.668137087446809 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    zbl = [0.0000000000000000 0.23677048111394497 0.54205398009286099 0.86821029420482398 1.1942305982289574 1.5021733114775246 1.7728696532854951 1.9811472656659215 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0735325177709347 2.0676124491406309 1.9538721645932462 1.7246580724060550 1.4199092524833361 1.0795744726067100 0.74157447412335120 0.43989215154552436 0.20287686084246737 5.1839272309835115E-002 1.2959818077458779E-002 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 6.4799090387293893E-003 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    sbl = [0.0000000000000000 0.25841299983870264 0.69309443916718450 1.3001178125821706 2.0803642302394660 3.0297111639838339 4.1405214334091553 5.4026510157169643 6.8043898816286141 8.3416218446335328 10.000552428355771 11.763006042166527 13.609672875455766 15.520320459958088 17.474015340908966 19.449352428355770 21.424689515802580 23.378384396753461 25.289031981255782 27.135698814545016 28.898152428355768 30.557093575238806 32.098527634941931 33.515875832618846 34.797462100476722 35.927263281032928 36.886781329295545 37.657175681161526 38.221334951814825 38.565663707653812 38.676246842534837 38.779972475720939 39.089407346005252 39.601363823293127 40.310232805285608 41.208247775075613 42.285569892732411 43.530395791690815 44.929086898905268 46.466318861910182 48.125249445632427 49.887703059443183 51.734369892732417 53.645017477234731 55.598712358185615 57.574049445632426 59.549386533079229 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    dybl = [0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    uinv = [0.0000000000000000 0.38031630119071819 0.65133000068646940 0.82986457636122890 0.94593900648830542 1.0273647497696450 1.0875710294281722 1.1379596292754699 1.1211770043036997 1.0410745224703368 1.0233748560204126 1.0153656859801528 1.0114146582262913 1.0093288558407103 1.0083771307464999 1.0082617907520570 1.0089589138885453 1.0107065580926284 1.0141479675815785 1.0212689130026895 1.0372367680217678 1.1144811750773000 1.1204432983136237 1.0698988951578923 0.99909801819802468 0.92090964221723726 0.84122659191418603 0.75823296164918463 0.65652375557614617 0.47983269527465694 0.48299181951594639 0.78041993927936115 0.87783297956297346 0.92544045471354175 0.95240359259102925 0.96862638059411565 0.97874472407836843 0.98523226140642761 0.98949362550703812 0.99235592689293306 0.99431919450472406 0.99569255108234722 0.99667107351754636 0.99738029809265483 0.99790253288612152 0.99829271154835397 0.99858812000216846 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]

    xl = zeros(nbldim)     # body x coordinates
    zl = zeros(nbldim)     # body z coordinates
    sl = zeros(nbldim)     # body + wake arc length  (BL coordinate)
    dyl = zeros(nbldim)     # body y-offset of edge-type tail
    ql = zeros(nbldim)     # inviscid edge velocity (in absence of delta*)

    nl, ilte = TASOPT.aerodynamics._axisymm_flow(xnose, xend, xblend1, xblend2,
        Sfuse, anose, btail, ifclose,
        Mach, nc, nbldim,
        xl, zl, sl, dyl, ql)

    @test nbl == nl
    @test iblte == ilte
    @test all(isapprox.(xl, xbl'))
    @test all(isapprox.(zl, zbl'))
    @test all(isapprox.(sl, sbl'))
    @test all(isapprox.(dyl, dybl'))
    @test all(isapprox.(ql, uinv'))

    # end _axisymm_flow

    #_axisymm_BL
    nbldim, nbl, iblte = [60 47 31]
    sbl = [0.0000000000000000 0.25841299983870264 0.69309443916718450 1.3001178125821706 2.0803642302394660 3.0297111639838339 4.1405214334091553 5.4026510157169643 6.8043898816286141 8.3416218446335328 10.000552428355771 11.763006042166527 13.609672875455766 15.520320459958088 17.474015340908966 19.449352428355770 21.424689515802580 23.378384396753461 25.289031981255782 27.135698814545016 28.898152428355768 30.557093575238806 32.098527634941931 33.515875832618846 34.797462100476722 35.927263281032928 36.886781329295545 37.657175681161526 38.221334951814825 38.565663707653812 38.676246842534837 38.779972475720939 39.089407346005252 39.601363823293127 40.310232805285608 41.208247775075613 42.285569892732411 43.530395791690815 44.929086898905268 46.466318861910182 48.125249445632427 49.887703059443183 51.734369892732417 53.645017477234731 55.598712358185615 57.574049445632426 59.549386533079229 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    bbl = [0.0000000000000000 1.4876728081089807 3.4058256034176799 5.4551261640898163 7.5035721481764730 9.4384332795128874 11.139268557067989 12.447915390991129 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 13.028389049617431 12.991192161382012 12.276540876679459 10.836346260450391 8.9215529527316466 6.7831664642886311 4.6594498399912689 2.7639239033344545 1.2747129112121085 0.32571575411203757 8.1428938528009392E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 4.0714469264004696E-002 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    rnbl = [0.25687167487068524 0.52549731623855278 0.77093224537019733 0.87388058682569725 0.92655898110869783 0.95779636867286155 0.97833423687623799 0.99298883438210350 0.99940498221588692 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 0.99999831036733711 0.99919998466767346 0.99280823486828462 0.97942256294624130 0.96239445494390818 0.94434885779005473 0.92736394348868623 0.91293247282059986 0.90204369867538925 0.92780848366116087 0.97883534072521383 0.99890404250382259 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 1.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    univ = [0.0000000000000000 0.38031630119071819 0.65133000068646940 0.82986457636122890 0.94593900648830542 1.0273647497696450 1.0875710294281722 1.1379596292754699 1.1211770043036997 1.0410745224703368 1.0233748560204126 1.0153656859801528 1.0114146582262913 1.0093288558407103 1.0083771307464999 1.0082617907520570 1.0089589138885453 1.0107065580926284 1.0141479675815785 1.0212689130026895 1.0372367680217678 1.1144811750773000 1.1204432983136237 1.0698988951578923 0.99909801819802468 0.92090964221723726 0.84122659191418603 0.75823296164918463 0.65652375557614617 0.47983269527465694 0.48299181951594639 0.78041993927936115 0.87783297956297346 0.92544045471354175 0.95240359259102925 0.96862638059411565 0.97874472407836843 0.98523226140642761 0.98949362550703812 0.99235592689293306 0.99431919450472406 0.99569255108234722 0.99667107351754636 0.99738029809265483 0.99790253288612152 0.99829271154835397 0.99858812000216846 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    Reunit = 7006067.9699355923
    Mach = 0.80000000000000004
    fex = 1.0300000000000000

    #Fortran results
    uebl = [0.0000000000000000 0.37728242338296181 0.64847633615516287 0.82695029534888154 0.94403149002962905 1.0255943149422753 1.0865499203351234 1.1365507804257704 1.1206431211650818 1.0440948790037765 1.0240821955777670 1.0158405647282294 1.0117347512309973 1.0095591700001416 1.0085471996187523 1.0083948355769166 1.0090823376403273 1.0108720300588245 1.0145032985652438 1.0221896394842940 1.0402657551827854 1.1059092729914484 1.1091568618457064 1.0543655732241883 0.97490990951276690 0.89368511298240294 0.83118999557894302 0.79419493273373520 0.77809105336189510 0.78142056972012308 0.79242533118539638 0.79772286109406265 0.81053258282277518 0.84123660173182713 0.88716224825623247 0.93244222532354892 0.96275025453959628 0.97871685796354857 0.98673875476417017 0.99107603733351624 0.99366250031367975 0.99532641373569586 0.99645251843379634 0.99724082638968936 0.99780393344085772 0.99820042263719722 0.99838556039722948 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    hkbl = [0.0000000000000000 1.3878006383330910 1.3464796562458770 1.2841329531938646 1.2957711548906594 1.2666086943831041 1.2667873219861920 1.2551499686354026 1.2703003810691702 1.2907324422017787 1.2593018975066568 1.2487806424396475 1.2423297359023246 1.2377736895702194 1.2341982436579075 1.2311759749694335 1.2284089245123413 1.2255808445615450 1.2221794255394689 1.2168451815800081 1.2059592263759851 1.1693954460989937 1.1814318659577450 1.2288851612693361 1.3020077541398134 1.4034862848531129 1.5185526081155398 1.6114897097201484 1.6546827444891374 1.6259707518242759 1.5786779780686202 1.5551481333557295 1.4996868800960224 1.4070761827361435 1.3136726651554176 1.2479933969997437 1.2114726616058633 1.1921668657173372 1.1806347627072433 1.1724280279402921 1.1658003248537754 1.1600708111266333 1.1549504781978848 1.1503025689488813 1.1460563390645413 1.1421769429235715 1.1386910005249200 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]
    phbl = [0.0000000000000000 3.0682260506125015E-005 4.7889944403344497E-004 2.4975159537457057E-003 7.6025183498762659E-003 1.7075491938605240E-002 3.1774804005064786E-002 5.1993342876357089E-002 7.5806565529789635E-002 9.9331571057468898E-002 0.12134723616883887 0.14314065217342706 0.16513947806483301 0.18732135100076602 0.20956926913881643 0.23173274901782004 0.25365105360216367 0.27516905047121493 0.29615615812688728 0.31655166909860810 0.33650147058369984 0.35720543526812448 0.37712570686668656 0.39265309496967099 0.40327292872315201 0.41058519440318958 0.41595841879182138 0.41988533035222408 0.42246231552153879 0.42384142776268363 0.42423434117475162 0.42468201186142845 0.42613631365539723 0.42790475253569571 0.42946904266257002 0.43066422096273016 0.43159729923635759 0.43240447120609804 0.43315867903042210 0.43388607014134134 0.43459219468942439 0.43527542653948964 0.43593220917558062 0.43655893064098156 0.43715256195251562 0.43771085141483057 0.43823251379636796 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000 0.0000000000000000]

    uei, dsi, thi, tsi,
    dci, cfi, cdi, cti,
    hki, phi = TASOPT.aerodynamics._axisymm_BL(nbldim, nbl, iblte,
        sbl, bbl, rnbl, univ,
        Reunit, Mach, fex)

    @test all(uei .≈ uebl')
    @test all(hki .≈ hkbl')
    @test all(phi .≈ phbl')

    #Fuse BL
    #TODO?

    
    #CMVf1 calculation
    ac = load_default_model(); 
    ac.fuselage.calculates_pitching_moment_volume = true
    size_aircraft!(ac, printiter=false)

    #ref value from TASOPT.jl, REPL execution
    CMVf1_test = 67.0315667
    @test ac.parg[igCMVf1] ≈ CMVf1_test
end

@testset "BL laminar dissipation" begin
  @testset "closure (Fortran DIL)" begin
    # Reference values evaluated from the Fortran DIL formulas (tasopt2.16 blsys.f:400)
    # in 40-digit decimal arithmetic, independently of the Julia code.
    # Columns: hk, Rθ, 2CD/H*, ∂/∂hk, ∂/∂Rθ. The first two rows are on the hk < 4
    # branch, the last two on the hk ≥ 4 branch, and hk = 4 is the switch point.
    ref = [
        2.5  500.0  4.52131666105232505e-4 -1.39816109052519186e-4  -9.04263332210465011e-7
        3.9 1000.0  2.07000006482669203e-4 -3.56546806183984770e-10 -2.07000006482669203e-7
        4.0  200.0  1.03500000000000000e-3  0.0                     -5.17500000000000000e-6
        5.0  300.0  6.84771241830065359e-4 -1.02524669998718442e-5  -2.28257080610021786e-6
        7.0 2000.0  9.73983050847457627e-5 -3.44728526285550129e-6  -4.86991525423728814e-8
    ]
    ld = TASOPT.aerodynamics.laminar_dissipation
    for (hk, rt, di_ref, di_hk_ref, di_rt_ref) in eachrow(ref)
        di, di_hk, di_rt = ld(hk, rt)
        @test di ≈ di_ref rtol = 1e-14
        @test di_hk ≈ di_hk_ref rtol = 1e-12 atol = 1e-20
        @test di_rt ≈ di_rt_ref rtol = 1e-14

        # partials against central finite differences (skipped at the hk = 4 kink,
        # where the slope is 0 but a central difference straddles both branches)
        ε = 1e-6
        fd_hk = (ld(hk + ε, rt)[1] - ld(hk - ε, rt)[1]) / (2ε)
        fd_rt = (ld(hk, rt + ε * rt)[1] - ld(hk, rt - ε * rt)[1]) / (2ε * rt)
        hk == 4.0 || @test di_hk ≈ fd_hk rtol = 1e-6 atol = 1e-13
        @test di_rt ≈ fd_rt rtol = 1e-6
    end

    # value and slope are continuous across the hk = 4 branch switch
    lo, hi = ld(4.0 - 1e-9, 300.0), ld(4.0 + 1e-9, 300.0)
    @test lo[1] ≈ hi[1] rtol = 1e-12
    @test lo[2] ≈ hi[2] atol = 1e-12
  end

  @testset "station vars, laminar branch" begin
    # Laminar branch of the station variables (is_laminar = true). Before the fix this
    # called an undefined function and threw UndefVarError.
    Reyn, Mach, fexcr = 7006067.9699355923, 0.8, 1.0
    x, θ, ue = 1.0, 0.0002, 0.95
    # (δs, hk, 2CD/H*) from the same independent evaluation of the Fortran blvar
    # chain; δs = 0.00052 gives hk < 4 and δs = 0.001 gives hk > 4.
    stations = [
        (0.00052, 2.28712293625379124, 1.81808278181889317e-4)  # Rθ = 1356.15,
        (0.001,   4.54177819543849179, 1.52293736895109506e-4)  # Rθ = 1356.15,
    ]
    vars(is_wake, θ, δs, ue; M = Mach) = TASOPT.aerodynamics._BL_station_vars(false, true,
        is_wake, Reyn, M, fexcr, x, θ, δs, ue)
    fd(f, v) = (f(v * (1 + 1e-6)) - f(v * (1 - 1e-6))) / (2e-6 * v)
    for (δs, hk_ref, di_ref) in stations
        out = @inferred vars(false, θ, δs, ue)
        hk = out[4]
        di, di_θ, di_δs, di_ue = out[20:23]
        @test hk ≈ hk_ref rtol = 1e-13
        @test di ≈ di_ref rtol = 1e-13

        # chain-rule partials against central finite differences of the station vars
        @test di_θ  ≈ fd(t -> vars(false, t, δs, ue)[20], θ)  rtol = 1e-5
        @test di_δs ≈ fd(d -> vars(false, θ, d, ue)[20], δs) rtol = 1e-5
        # The ue partial at M = 0.8 needs the exact msq_ue (the Fortran blvar drops a
        # 1/trat factor, which is exact only at M = 0), and is also checked at M = 0.
        @test di_ue ≈ fd(u -> vars(false, θ, δs, u)[20], ue) rtol = 1e-5
        di_ue0 = vars(false, θ, δs, ue; M = 0.0)[23]
        @test di_ue0 ≈ fd(u -> vars(false, θ, δs, u; M = 0.0)[20], ue) rtol = 1e-5

        # laminar wake: the Fortran blvar doubles the dissipation (wfac = 2)
        wake = vars(true, θ, δs, ue)
        @test wake[20] ≈ 2di rtol = 1e-14
        @test wake[16] == 0.0   # cf = 0 in the wake
    end
  end
end

@testset "BL closures and station partials" begin
  # The BL closures were ported from Fortran (tasopt2.16 blsys.f) with hand-coded
  # partials. Check every partial against ForwardDiff of the value, which also catches
  # slips that are in the Fortran too.
  AE = TASOPT.aerodynamics
  # Each closure f(p...) returns (value, ∂value/∂p[1], ∂value/∂p[2], ...).
  function check_partials(f, p)
    out = f(p...)
    for k in eachindex(p)
      ad = ForwardDiff.derivative(t -> f(Base.setindex(p, t, k)...)[1], p[k])
      @test out[k+1] ≈ ad rtol = 1e-12
    end
  end

  @testset "closures" begin
    for H in (1.3, 2.0, 3.5), msq in (0.0, 0.64)
      check_partials(AE.hkin, (H, msq))
    end
    for hk in (2.5, 3.9, 5.0, 7.0), rt in (300.0, 2000.0)
      check_partials(AE.laminar_dissipation, (hk, rt))
      check_partials(AE.DILW, (hk, rt))
    end
    for hk in (2.0, 4.0, 5.0, 7.0)
      check_partials(AE.hsl, (hk, 500.0, 0.3))
    end
    for hk in (2.0, 5.0, 6.0, 8.0), rt in (300.0, 2000.0)
      check_partials(AE.cfl, (hk, rt, 0.3))
    end
    check_partials(AE.DIT, (1.6, 0.3, 0.003, 0.02))
    for hk in (1.3, 2.5), msq in (0.0, 0.5)
      check_partials(AE.hct, (hk, msq))
    end
    for hk in (1.2, 1.5, 2.5, 4.0), rt in (300.0, 5000.0, 1.0e5), msq in (0.0, 0.5)
      check_partials(AE.cft, (hk, rt, msq))
    end
    # hst: Rθ on both sides of the Rθ = 200 and Rθ = 400 switches, Hk on the attached
    # (Hk < H₀) and separated (Hk ≥ H₀) branches.
    @testset "hst Rθ = $rt" for rt in (150.0, 300.0, 2000.0, 1.0e5)
      H₀ = rt > 400 ? 3 + 400 / rt : 4.0
      for hk in (1.3, 2.5, H₀ + 0.3, H₀ + 2.0), msq in (0.0, 0.5)
        check_partials(AE.hst, (hk, rt, msq))
      end
    end
  end

  # _BL_station_vars returns (h, h_θ, h_δs, then (q, q_θ, q_δs, q_ue) for
  # q = hk, hc, hs, cf, 𝒟).
  vidx = (1, 4, 8, 12, 16, 20)
  Reyn, fexcr, x = 7.0e6, 1.03, 10.0
  # (θ, δs): attached (H = 1.3, 2.5) and separated (H = 6, Hk > H₀)
  states = ((2e-3, 2.6e-3), (2e-3, 5e-3), (2e-4, 1.2e-3))

  @testset "station vars, Mach = $Mach" for Mach in (0.0, 0.8)
    for is_laminar in (false, true), is_wake in (false, true), (θ, δs) in states
      ue = 0.95
      f(t, d, u) = AE._BL_station_vars(false, is_laminar, is_wake, Reyn, Mach, fexcr, x, t, d, u)
      out = f(θ, δs, ue)
      for vi in vidx
        @test out[vi+1] ≈ ForwardDiff.derivative(t -> f(t, δs, ue)[vi], θ) rtol = 1e-12
        @test out[vi+2] ≈ ForwardDiff.derivative(d -> f(θ, d, ue)[vi], δs) rtol = 1e-12
        vi == 1 && continue   # h = δs/θ does not depend on ue
        @test out[vi+3] ≈ ForwardDiff.derivative(u -> f(θ, δs, u)[vi], ue) rtol = 1e-12 atol = 1e-14
      end
    end
  end

  # _BL_station_system: aa = ∂rr/∂(θ, δs, ue) and bb = ∂rr/∂(θ, δs, ue) of the previous
  # station, with the station variables recomputed from (θ, δs, ue). The system fills
  # Float64 static arrays, so the reference is a central finite difference.
  function station_residual(simi, wake, direct, Mach, s, sm)
    v  = AE._BL_station_vars(simi, false, wake, Reyn, Mach, fexcr, x, s...)
    vm = AE._BL_station_vars(false, false, wake, Reyn, Mach, fexcr, 0.9x, sm...)
    AE._BL_station_system(simi, false, wake, direct, Mach, 1.0, 2.9,
                          x, 12.0, 0.9, s..., v...,
                          0.9x, 12.5, 0.95, sm..., vm...)
  end
  @testset "station system, Mach = $Mach" for Mach in (0.0, 0.8)
    s, sm = (2.1e-3, 2.9e-3, 0.97), (2.0e-3, 2.7e-3, 0.99)
    for simi in (false, true), wake in (false, true), direct in (true, false)
      aa, bb, rr = station_residual(simi, wake, direct, Mach, s, sm)
      for j in 1:3
        step(p, δ) = Base.setindex(p, p[j] * (1 + δ), j)
        h = 1e-6
        fd_a = (station_residual(simi, wake, direct, Mach, step(s, h), sm)[3] -
                station_residual(simi, wake, direct, Mach, step(s, -h), sm)[3]) / (2h * s[j])
        fd_b = (station_residual(simi, wake, direct, Mach, s, step(sm, h))[3] -
                station_residual(simi, wake, direct, Mach, s, step(sm, -h))[3]) / (2h * sm[j])
        for k in 1:3
          @test aa[k, j] ≈ fd_a[k] rtol = 1e-6 atol = 1e-9 * maximum(abs, aa[k, :])
          # the self-similar station does not depend on the previous one
          simi || @test bb[k, j] ≈ fd_b[k] rtol = 1e-6 atol = 1e-9 * maximum(abs, bb[k, :])
        end
      end
    end
  end
end

@testset "BL Newton system: block elimination" begin
  # The global viscous/inviscid Newton pass of _axisymm_BL solves a 3n×3n system with
  # unknowns (δθᵢ, δmᵢ, δueᵢ) per station: rows 1–2 of station i couple stations i-1 and
  # i only, row 3 (ue - uvis) is dense in the m columns, and station 1 is frozen
  # (identity rows). _BL_block_solve! uses that structure; check it against dense LU.
  AE = TASOPT.aerodynamics

  # deterministic pseudo-random numbers in [-1, 1]
  prand(k) = sin(12.9898 * k + 78.233 * sin(0.37 * k))
  # A matrix with the _axisymm_BL Newton-system structure. `singular_at = s` makes the
  # 2×2 (θ, ue) block of station s exactly singular.
  function structured_system(n; singular_at = 0)
    N = 3n
    A, r = zeros(N, N), zeros(N)
    A[1, 1] = A[2, 2] = A[3, 3] = 1.0
    r[1:3] .= (0.1, -0.2, 0.05)       # generic: station 1 need not be zero
    c = 0
    for i in 2:n
      k = 3(i - 1)
      for row in 1:2, col in 1:3
        A[k+row, k+col]   = prand(c += 1) + (row == col ? 3.0 : 0.0) + (row == 2 && col == 3 ? 3.0 : 0.0)
        A[k+row, k+col-3] = 0.5 * prand(c += 1)
      end
      if i == singular_at                # (θ, ue) columns of rows 1–2 linearly dependent
        A[k+2, k+1] = 2A[k+1, k+1]
        A[k+2, k+3] = 2A[k+1, k+3]
      end
      for j in 1:n
        A[k+3, 3(j-1)+2] = -0.3 * prand(c += 1) / (1 + abs(i - j))
      end
      A[k+3, k+3] = 1.0
      r[k+1:k+3] .= (prand(c += 1), prand(c += 1), prand(c += 1))
    end
    return A, r
  end

  @testset "structured system, n = $n" for n in (2, 3, 10, 47, 60)
    A, r = structured_system(n)
    ws = AE._BLBlockWork(n)
    z = copy(r)
    @test AE._BL_block_solve!(z, copy(A), n, ws)    # true: block elimination was used
    zd = A \ r
    @test maximum(abs, z - zd) <= 1e-12 * maximum(abs, zd)
    @test maximum(abs, A * z - r) <= 1e-12 * maximum(abs, r)
  end

  @testset "near-singular 2×2 block falls back to dense LU" begin
    n = 12
    A, r = structured_system(n; singular_at = 5)
    z = copy(r)
    @test !AE._BL_block_solve!(z, copy(A), n, AE._BLBlockWork(n))   # false: dense LU used
    @test maximum(abs, A * z - r) <= 1e-10 * maximum(abs, r)
  end

  # _axisymm_BL inputs built from fuselage geometry the way fuselage_drag! builds them
  function fuse_bl_inputs(; Mach = 0.8, nc = 30, edge = false, xend = 37.7952,
                          xblend1 = 6.096, xblend2 = 29.5656, btail = 2.0)
    nbldim = 60
    xl, zl, sl, dyl, ql = (zeros(nbldim) for _ in 1:5)
    nbl, ite = AE._axisymm_flow(0.0, xend, xblend1, xblend2, 13.507394174276255,
                                1.65, btail, Int(edge), Mach, nc, nbldim, xl, zl, sl, dyl, ql)
    bl = edge ? 2π .* zl .+ 4 .* dyl : 2π .* zl
    rn = zeros(nbldim)
    for i in 2:nbl-1
      dxm, dzm, dsm = xl[i] - xl[i-1], zl[i] - zl[i-1], sl[i] - sl[i-1]
      dxp, dzp, dsp = xl[i+1] - xl[i], zl[i+1] - zl[i], sl[i+1] - sl[i]
      dxo = (dxm / dsm) * dsp + (dxp / dsp) * dsm
      dzo = (dzm / dsm) * dsp + (dzp / dsp) * dsm
      rn[i] = dxo / sqrt(dxo^2 + dzo^2)
    end
    rn[1] = max(rn[2] - (0.5(sl[3] + sl[2]) - sl[1]) * (rn[3] - rn[2]) / (sl[3] - sl[2]), 0.0)
    rn[nbl] = max(rn[nbl-1] + (sl[nbl] - 0.5(sl[nbl-1] + sl[nbl-2])) *
                  (rn[nbl-1] - rn[nbl-2]) / (sl[nbl-1] - sl[nbl-2]), 0.0)
    return nbldim, nbl, ite, sl, bl, rn, ql
  end

  cases = (
    ("default, M = 0.8",      (;)),
    ("M = 0",                 (; Mach = 0.0)),
    ("M = 0.6",               (; Mach = 0.6)),
    ("M = 0.85",              (; Mach = 0.85)),
    ("tapers to edge",        (; edge = true)),
    ("n = 26",                (; nc = 16)),
    ("n = 59",                (; nc = 38)),
    ("short tail, Hk ≈ 2.4",  (; xblend2 = 34.0)),
    ("short tail, Hk > 2.9",  (; xblend2 = 35.0)),   # inverse mode in the initial march
  )
  @testset "_axisymm_BL, $label" for (label, kw) in cases
    Mach = get(kw, :Mach, 0.8)
    nbldim, nbl, ite, sl, bl, rn, ql = fuse_bl_inputs(; kw...)
    args = (nbldim, nbl, ite, sl, bl, rn, ql, 7006067.9699355923, Mach, 1.03)
    blk = AE._axisymm_BL(args...)
    dns = AE._axisymm_BL(args...; dense_solve = true)
    # uₑ, δ*, θ, θ*, δ**, cf, c𝒟, cτ, Hk, Φ: same to rounding
    for (b, d) in zip(blk, dns)
      @test maximum(abs, b - d) <= 1e-11 * maximum(abs, d)
    end
    @test all(isfinite, blk[1])
    nbl > ite && @test all(iszero, blk[6][ite+1:nbl])   # cf = 0 in the wake
  end
end

@testset "aeroperf sweep" begin
    #loading known results, hardcoded; can be re-generated via aeroperf_sweep run in REPL
    test_results = (
        CLs = [0.0, 0.8],
        CDs = [0.13614022933824074, 0.07096802425920105],
        LDs = [0.0, 11.272682427766467],
        CLhs = [-0.03504971960044751, 0.031985774363004225],
        CDis = [0.0, 0.022084821719897365],
        CDwings = [0.12319577583754379, 0.03593874903860673],
        CDfuses = [0.006779753907171779, 0.006779753907171779],
        CDhtails = [0.0025503687249054377, 0.0025503687249054377],
        CDvtails = [0.0017228179947728997, 0.0017228179947728997],
        CDothers = [0.001891512873846836, 0.001891512873846836],
        clpos = [0.013754372953093726, 0.8876652409539934],
        clpss = [0.01716252120582795, 1.10761672476844],
        clpts = [0.012713905510437195, 0.8205166485506765],
        cdfss = [0.004870853045517493, 0.0042903934423377515],
        cdpss = [0.1489760335583573, 0.06713785293651943],
        cdwss = [0.0, 0.0],
        cdss = [0.1538468866038748, 0.07142824637885718]
    )

    #get basic execution of default model
    ac1 = load_default_model()
    size_aircraft!(ac1, printiter=false)
    results = aeroperf_sweep(ac1, [0.0, 0.8])

    #compare the two for approximate equality
    for field in keys(test_results)
        @test isapprox(results[field], test_results[field]; atol=1e-5)
    end
end
