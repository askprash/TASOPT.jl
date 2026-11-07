#Sample parameters
ac = load_default_model()

ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
ac.options.ifuel = ifuel
ac.para[iaalt, ipcruise1,1] = z
ac.para[iaMach, ipcruise1,1] = Mair
ac.parg[igxftank] = xftank
ac.fuse_tank.tank_count = 1
ac.parg[igWfuel] = 1e5

fuse_tank = ac.fuse_tank
fuse_tank.fueltype = "LH2"
fuse_tank.placement = "front"

fuse_tank.clearance_fuse = 0.1
TSL = 288.2
fuse_tank.TSLtank = [TSL]

fuse_tank.pvent = 2e5
fuse_tank.ARtank = 2.0
fuse_tank.inner_material = TASOPT.StructuralAlloy("Al-2219-T87")

fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = [ThermalInsulator("rohacell41s"), ThermalInsulator("rohacell41s"), ThermalInsulator("polyurethane27")]
fuse_tank.sizes_insulation = true
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.boiloff_rate = 0.15

fuse_tank.ullage_frac = 0.1
fuse_tank.qfac = 1.3
fuse_tank.theta_inner = 1.0
fuse_tank.ftankadd = 0.1
fuse_tank.ew = 0.9

fuse_tank.outer_material = StructuralAlloy("Al-2219-T87")
fuse_tank.theta_outer = [1.0, 2.0]
fuse_tank.Ninterm = 1.0

#Calculate fuel temperature and density as a function of pressure
β0 = 1 - fuse_tank.ullage_frac
fuel_mix = TASOPT.SaturatedMixture(fuse_tank.fueltype, fuse_tank.pvent, β0)

Tsat = fuel_mix.liquid.Tsat
ρl = fuel_mix.liquid.ρ
ρg = fuel_mix.gas.ρ
hvap = fuel_mix.hvap

fuse_tank.rhofuel = ρl
fuse_tank.rhofuelgas = ρg
fuse_tank.Tfuel = Tsat
fuse_tank.hvap = hvap

#Fuselage
fuse = ac.fuselage
fuse.layout.cross_section.radius = 2.5
fuse.layout.cross_section.bubble_lower_downward_shift = 0.3

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin

        TASOPT.tanksize!(ac, 1)
        inner_tank_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        @test inner_tank_mech.Wtank       ≈ 60109.124334566346
        @test inner_tank_mech.Winsul_sum  ≈ 38147.81559388217
        @test inner_tank_mech.Vfuel       ≈ 166.77327116515787
        @test inner_tank_mech.Shead_insul ≈ [16.42133035129676, 20.311090752967278, 24.620380536430755, 29.34797626786876]
        @test inner_tank_mech.Rtank_outer ≈ 1.8743852420176998
        @test inner_tank_mech.l_tank      ≈ 15.771807531012154
        @test inner_tank_mech.l_cyl       ≈ 12.842701653674803
        
        outputs_thermal = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, TSL, Mair, xftank, ifuel)

        outputs_thermal_check = (1835.8608916011449,)

        for i in 1:length(outputs_thermal)
            @test outputs_thermal[i] ≈ outputs_thermal_check[i]
        end
    end
    
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.sizes_insulation = false

    @testset "Vacuum insulation" begin
        TASOPT.tanksize!(ac, 1)

        inner_tank_vac = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        @test inner_tank_vac.Wtank       ≈ 24229.845334360776
        @test inner_tank_vac.Winsul_sum  ≈ 0.0
        @test inner_tank_vac.Vfuel       ≈ 166.77327116515787
        @test inner_tank_vac.Shead_insul ≈ [25.483812169139977, 27.21328988831087]
        @test inner_tank_vac.Rtank_outer ≈ 2.335
        @test inner_tank_vac.l_tank      ≈ 9.994915110929762
        @test inner_tank_vac.l_cyl       ≈ 7.525566077704729

        Winnertank = inner_tank_vac.Wtank
        l_cyl      = inner_tank_vac.l_cyl
        Winner_tot = Winnertank + fuse_tank.Wfuelintank

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.CryoTank.optimize_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl)
        Ninterm_check = 14.025390625
        @test Ninterm ≈ Ninterm_check

        outer_tank_vac = TASOPT.CryoTank.size_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl, Ninterm_check)
        @test outer_tank_vac.Wtank   ≈ 87481.40843690369
        @test outer_tank_vac.Wcyl    ≈ 31381.827647692957
        @test outer_tank_vac.Whead   ≈ 13786.613724211644
        @test outer_tank_vac.Wstiff  ≈ 20573.498028341648
        @test outer_tank_vac.Souter  ≈ 171.6223477616751
        @test outer_tank_vac.Shead   ≈ 26.90237940128641
        @test outer_tank_vac.Scyl    ≈ 117.81758895910225
        @test outer_tank_vac.t_cyl   ≈ 0.009560503271690743
        @test outer_tank_vac.t_head  ≈ 0.018394143361195915
        @test outer_tank_vac.l_outer ≈ 9.96235436442712

    end
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.CryoTank.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1749999999999998, 0.09494892091249901]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM[i] ≈ outputs_bendM_check[i]
        end

        θ1 = 45 * pi /180
        θ2 = 135 * pi /180
        outputs_bendM_outer = TASOPT.CryoTank.stiffeners_bendingM_outer(θ1, θ2)
        outputs_bendM_outer_check = [pi/2, 0.6166583566771067]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM_outer[i] ≈ outputs_bendM_outer_check[i]
        end

        perim_vessel = fuse.layout.cross_section.perimeter

        Wstiff = TASOPT.CryoTank.stiffener_weight("outer", 7e4, fuse.layout.radius, perim_vessel, fuse_tank.outer_material.UTS / 4, 
        fuse_tank.outer_material.ρ, θ1, θ2, 10.0, 5.0, fuse_tank.outer_material.E)
        Wstiff_check = 1717.0752091513864
        
        @test Wstiff ≈ Wstiff_check
    end

    @testset "Thermal models" begin
        outputs_h = TASOPT.CryoTank.tank_heat_coeff(21.0, ifuel, fuse_tank.Tfuel, 5.0)
        outputs_h_check = (133.43555789562762)
        for i in 1:length(outputs_h)
            @test outputs_h[i] ≈ outputs_h_check[i]
        end

        outputs_free = TASOPT.CryoTank.freestream_heat_coeff(z, TSL, Mair, xftank, 240.0)
        outputs_free_check = (91.30101055991513, 218.06145060705106, 243.26644796151263)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        outputs_natural = TASOPT.CryoTank.freestream_heat_coeff(0.0, TSL, 0.0, xftank, 287.0, fuse.layout.radius)
        outputs_natural_check = (1.6795624846472137, 288.2, 288.2)
        for i in 1:length(outputs_natural)
            @test outputs_natural[i] ≈ outputs_natural_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.3565313051181731

        @test Rvac ≈ Rvac_check
    end
end

@testset "Tank thermal solves vs NLsolve reference" begin
    # Fresh copy of the LH2 tank set up at the top of this file: tanksize! stores the tank
    # position in parg, which feeds the next insulation sizing, so each case starts anew.
    function lh2_tank_aircraft(t_insul, materials)
        ac = load_default_model()
        ac.options.ifuel = 40
        ac.para[iaalt, ipcruise1, 1] = 11e3
        ac.para[iaMach, ipcruise1, 1] = 0.8
        ac.parg[igxftank] = 15.0
        ac.parg[igWfuel] = 1e5
        ft = ac.fuse_tank
        ft.tank_count = 1
        ft.fueltype = "LH2"
        ft.placement = "front"
        ft.clearance_fuse = 0.1
        ft.TSLtank = [288.2]
        ft.pvent = 2e5
        ft.ARtank = 2.0
        ft.inner_material = StructuralAlloy("Al-2219-T87")
        ft.ullage_frac = 0.1
        ft.qfac = 1.3
        ft.theta_inner = 1.0
        ft.ftankadd = 0.1
        ft.ew = 0.9
        ft.outer_material = StructuralAlloy("Al-2219-T87")
        ft.theta_outer = [1.0, 2.0]
        ft.Ninterm = 1.0
        mix = TASOPT.SaturatedMixture(ft.fueltype, ft.pvent, 1 - ft.ullage_frac)
        ft.rhofuel = mix.liquid.ρ
        ft.rhofuelgas = mix.gas.ρ
        ft.Tfuel = mix.liquid.Tsat
        ft.hvap = mix.hvap
        ac.fuselage.layout.cross_section.radius = 2.5
        ac.fuselage.layout.cross_section.bubble_lower_downward_shift = 0.3
        ft.t_insul = copy(t_insul)
        ft.material_insul = ThermalInsulator.(materials)
        return ac
    end

    # The references are NLsolve (trust region) solutions at ftol = 1e-10, xtol = 0.
    # NLsolve at the production tolerances differs from them by at most 3e-10 relative
    # in Q and 1e-12 relative in the insulation thickness.
    @testset "tankWthermal" begin
        insulations = Dict("foam3" => ([0.1, 0.1, 0.1], ["rohacell41s", "rohacell41s", "polyurethane27"]),
                           "foam1" => ([0.12], ["polyurethane32"]),
                           "vacuum" => ([0.065], ["vacuum"]))
        # (insulation, z, Mair, TSL, Q); Mair = 0 uses the natural-convection model
        cases = [
            ("foam3", 0.0, 0.00, 288.2, 3533.3349591636729),
            ("foam3", 0.0, 0.25, 288.2, 3683.6452638964229),
            ("foam3", 0.0, 0.25, 308.2, 4076.2266233436658),
            ("foam3", 5000.0, 0.60, 288.2, 3326.3362909575562),
            ("foam3", 11000.0, 0.80, 288.2, 2787.7450389284791),
            ("foam3", 12500.0, 0.85, 288.2, 2818.6325215139095),
            ("foam1", 0.0, 0.00, 288.2, 9226.8934410410875),
            ("foam1", 0.0, 0.25, 288.2, 10191.752554412769),
            ("foam1", 0.0, 0.25, 308.2, 11586.297984103241),
            ("foam1", 5000.0, 0.60, 288.2, 8915.1966333320379),
            ("foam1", 11000.0, 0.80, 288.2, 7108.4523931699032),
            ("foam1", 12500.0, 0.85, 288.2, 7204.0402924024602),
            ("vacuum", 0.0, 0.00, 288.2, 2262.4835264279382),
            ("vacuum", 0.0, 0.25, 288.2, 2448.2305681526764),
            ("vacuum", 0.0, 0.25, 308.2, 3034.9597269352598),
            ("vacuum", 5000.0, 0.60, 288.2, 1987.9939529489443),
            ("vacuum", 11000.0, 0.80, 288.2, 1444.6108193309101),
            ("vacuum", 12500.0, 0.85, 288.2, 1471.489265787596),
        ]
        for (name, (t_insul, materials)) in insulations
            ac = lh2_tank_aircraft(t_insul, materials)
            ac.fuse_tank.sizes_insulation = false
            TASOPT.tanksize!(ac, 1) # tank geometry for this insulation
            for (case, z, Mair, TSL, Q) in cases
                case == name || continue
                @test TASOPT.CryoTank.tankWthermal(ac.fuselage, ac.fuse_tank, z, TSL, Mair, 15.0, 40) ≈ Q rtol = 1e-9
            end
        end
    end

    @testset "tanksize! insulation sizing" begin
        # (initial thicknesses, materials, layers sized, boil-off rate in %/h, sized thicknesses)
        cases = [
            ([0.05, 0.05, 0.05], ["polyurethane32", "polyurethane32", "polyurethane32"], [1, 2, 3], 0.4, [0.05831367625836618, 0.05831367625836618, 0.05831367625836618]),
            ([0.15, 0.15, 0.15], ["rohacell41s", "rohacell41s", "polyurethane27"], [1, 2, 3], 0.15, [0.17520491932727467, 0.17520491932727467, 0.17520491932727467]),
            ([0.15, 0.15, 0.15], ["rohacell41s", "rohacell41s", "polyurethane27"], [2], 0.2, [0.15, 0.06034262289978137, 0.15]),
            ([0.15, 0.15, 0.15], ["rohacell41s", "rohacell41s", "polyurethane27"], [3], 0.15, [0.15, 0.15, 0.2421410067143236]),
            ([0.05], ["polyurethane32"], [1], 0.8, [0.08542341281734236]),
        ]
        for (t_insul, materials, iinsuldes, boiloff, t_sized) in cases
            ac = lh2_tank_aircraft(t_insul, materials)
            ac.fuse_tank.sizes_insulation = true
            ac.fuse_tank.iinsuldes = iinsuldes
            ac.fuse_tank.boiloff_rate = boiloff
            TASOPT.tanksize!(ac, 1)
            for (t, tref) in zip(ac.fuse_tank.t_insul, t_sized)
                @test t ≈ tref rtol = 1e-9
            end
        end
    end
end
