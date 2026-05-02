"""
    function tfoper!(gee, M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, eng_has_BLI_cores,
      pid, pib, pifn, pitn,
      Gearf,
      pi_fan_des, pi_lpc_des, pi_hpc_des, pi_hpt_des, pi_lpt_des,
      mb_fan_des, mb_lpc_des, mb_hpc_des, mb_hpt_des, mb_lpt_des,
      Nb_fan_des, Nb_lpc_des, Nb_hpc_des, Nb_hpt_des, Nb_lpt_des,
      A2, A2_5, A5, A7,
      opt_calc_call,
      Ttf, ifuel, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
      mofft, Pofft,
      Tt25off, pt25off,
      epsl, epsh,
      opt_cooling,
      Mtexit, dTstrk, StA, efilm, tfilm,
      M4a, ruc,
      ncrowx, ncrow,
      epsrow, Tmrow,
      Feng,
      M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt8, mcore, M2_5)

Turbofan operation routine. Calculation procedure follows that of Kerrebrock, but
the usual gas property formulas are replaced by function calls, which can therefore
implement more general gas models. In addition, a turbine cooling model is added.

The gas routines are described in [Gas Calculations](@ref)

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `Tref`:    reference temperature for corrected mass flow and speed
    - `pref`:    reference pressure for corrected mass flow
    - `Phiinl`:  inlet ingested dissipation  `Phi_inl`
    - `eng_has_BLI_cores`:
      - `false`: core in clear flow
      - `true`: core sees `Phiinl`
    - `pid`:     diffuser pressure ratio  ( = `pt2/pt0`)
    - `pib`:     burner   pressure ratio  ( = `pt4/pt3`)
    - `pifn`:    fan     nozzle pressure ratio  ( = `pt18/pt9.9`)
    - `pitn`:    turbine nozzle pressure ratio  ( = `pt8/pt4.9`)
    - `Gearf`:   fan gear ratio  ( = Nl/Nf )
    - `pi_fan_des`:    design fan pressure ratio  ( = `pt13/pt2`)
    - `pi_lpc_des`:   design LPC pressure ratio  ( = `pt2_5/pt2a`)
    - `pi_hpc_des`:   design HPC pressure ratio  ( = `pt3/pt2_5`)
    - `pi_hpt_des`:   design HPT pressure ratio  ( = `pt4_5/pt4_1`)
    - `pi_lpt_des`:   design LPT pressure ratio  ( = `pt5/pt4_5`)
    - `mb_fan_des`:    design corrected fan mass flow ( = `mf*sqrt(Tt2/Tref)/(pt2/pref)` )
      where `mf = mc*BPR`
    - `mb_lpc_des`:   design corrected LPC mass flow ( = `mc*sqrt(Tt2a/Tref)/(pt2a/pref)` )
    - `mb_hpc_des`:   design corrected HLC mass flow ( = `mc*sqrt(Tt2_5/Tref)/(pt2_5/pref)` )
    - `mb_hpt_des`:   design corrected HPT mass flow ( = `mt*sqrt(Tt4_1/Tref)/(pt4_1/pref)` )
      where `mt = mc*(1+ff)`
    - `mb_lpt_des`:   design corrected LPT mass flow ( = `mt*sqrt(Tt4_5/Tref)/(pt4_5/pref)` )
    - `Nb_fan_des`:    design corrected fan speed ( = `Nf/sqrt(Tt2/Tref)` )
    - `Nb_lpc_des`:   design corrected LPC speed ( = `Nl/sqrt(Tt2a/Tref)` )
    - `Nb_hpc_des`:   design corrected HPC speed ( = `Nh/sqrt(Tt2_5/Tref)` )
    - `Nb_hpt_des`:   design corrected HPT speed ( = `Nh/sqrt(Tt4_1/Tref)` )
    - `Nb_lpt_des`:   design corrected LPT speed ( = `Nl/sqrt(Tt4_5/Tref)` )
    - `A2`:      fan-face area [m^2]
    - `A2_5`:     HPC-face area [m^2]
    - `A5`:      core nozzle area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `opt_calc_call::CalcMode.T`:
      - `CalcMode.FixedTt4OffDes`: `Tt4` is specified
      - `CalcMode.FixedFeOffDes`: `Feng` is specified
    - `Tt4`:     turbine-inlet total temperature [K]
    - `Ttf`:     fuel temperature entering combustor
    - `ifuel`:   fuel index, see function [`gasfun`](@ref)
    - `hvap`:    fuel enthalpy of vaporization (J/kg)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    max fan polytropic efficiency
    - `eplc0`:   LPC max polytropic efficiency
    - `ephc0`:   HPC max polytropic efficiency
    - `epht0`:   HPT max polytropic efficiency
    - `eplt0`:   LPT max polytropic efficiency

    - `mofft`:    mass flow offtake at LPC discharge station 2.5
    - `Pofft`:    low spool power offtake
    - `Tt25off`:     offtake air discharge total temperature
    - `pt25off`:     offtake air discharge total pressure
    - `epsl`:    low  spool power loss fraction
    - `epsh`:    high spool power loss fraction

    - `opt_cooling::CoolingOpt.T`: turbine cooling model
      - `CoolingOpt.NoCooling`: no cooling, ignore all cooling parameters below
      - `CoolingOpt.FixedCoolingFlowRatio`: cooling flow ratios `epsrow` are inputs; compute `Tmrow`
      - `CoolingOpt.FixedTmetal`: metal temperatures `Tmrow` are inputs; compute `epsrow`
    - `Mtexit`:   turbine blade-row exit Mach, for setting temperature drops
    - `Tmetal`:   specified metal temperature  [K], used only if `opt_cooling=CoolingOpt.FixedTmetal`
    - `dTstrk`:   hot-streak temperature delta {K}, used only if `opt_cooling=CoolingOpt.FixedTmetal`
    - `StA`:      area-weighted Stanton number    , used only if `opt_cooling=CoolingOpt.FixedTmetal`
    - `M4a`:      effective Mach at cooling-flow outlet (start of mixing)
    - `ruc`:      cooling-flow outlet velocity ratio, `u/ue`
    - `ncrowx`:      dimension of `epsrow` array
    - `ncrow`:       number of blade rows requiring cooling
    - `epsrow(.)`: input specified cooling-flow bypass ratio if `opt_cooling=CoolingOpt.FixedCoolingFlowRatio`;
      output resulting cooling-flow bypass ratio if `opt_cooling=CoolingOpt.FixedTmetal`.
    - `Tmrow(.)`: input specified metal temperature [K] if `opt_cooling=CoolingOpt.FixedTmetal`;
      output resulting metal temperature [K] if `opt_cooling=CoolingOpt.FixedCoolingFlowRatio`

    **Output:**
    - `epsrow(.)`:   see above
    - `Tmrow(.)`:    see above
    - `TSFC`:    thrust specific fuel consumption = `mdot_fuel g / F`   [1/s]
    - `Fsp`:     specific thrust  = `F / (mdot u0) = F / ((1+BPR) mdot_core u0)`
    - `hfuel`:   fuel heating value   [J / kg K]
    - `ff`:      fuel mass flow fraction  =  `mdot_fuel / mdot_core`
    - `Feng`:    net effective thrust  = `(PK_inl+PK_out-Phi_jet)/u0  =  sum(mdot u)`
    - `mcore`:   core mass flow = `mdot_core`  [kg/s]
    - `BPR`:     bypass ratio   = `mdot_fan/mdot_core`
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= `dh/dT`)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:      static temperature
    - `u?`:      velocity
    - `etaf`:    fan          overall efficiency
    - `etac`:    compressor   overall efficiency
    - `etatf`:   fan-turbine  overall efficiency
    - `etatc`:   comp-turbine overall efficiency
    - `Lconv`:   `true` if convergence was successful, `false` otherwise

    The "?" symbol denotes the station index:
    - 0:    freestream
    - 1.8:  fan face outside of casing BLs
    - 1.9:  fan face over LPC portion
    - 2:    fan face over fan portion
    - 2.1:  fan exit, precooler inlet
    - 1.9c: precooler outlet, LPC inlet
    - 2.5:  LPC exit, intercooler inlet
    - 2.5c: intercooler exit, HPC inlet
    - 3:    compressor exit
    - 4:    combustor exit before cooling air addition
    - 4.1:  turbine inlet after cooling air addition
    - 4.5:  HPT exit, LPT inlet
    - 4.9:  LPT exit, regenerative cooler inlet
    - 4.9c: regenerative cooler outlet
    - 5:    core nozzle
    - 6:    core flow downstream
    - 7:    fan nozzle
    - 8:    fan flow downstream
"""
function tfoper!(gee, M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, eng_has_BLI_cores,
      pid, pib, pifn, pitn,
      Gearf,
      pi_fan_des, pi_lpc_des, pi_hpc_des, pi_hpt_des, pi_lpt_des,
      mb_fan_des, mb_lpc_des, mb_hpc_des, mb_hpt_des, mb_lpt_des,
      Nb_fan_des, Nb_lpc_des, Nb_hpc_des, Nb_hpt_des, Nb_lpt_des,
      A2, A2_5, A5, A7,
      opt_calc_call,
      Ttf, ifuel, hvap, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
      mofft, Pofft,
      Tt25off, pt25off,
      epsl, epsh,
      opt_cooling,
      Mtexit, dTstrk, StA, efilm, tfilm,
      fc0, epht_fc,
      M4a, ruc,
      ncrowx, ncrow,
      epsrow, Tmrow,
      Feng,
      M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt8, mcore, M2_5, 
      Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
      Δp_PreC, Δp_InterC, Δp_Regen)

      #---- ncrowy must be at least as big as ncrowx defined in index.inc
      ncrowy = 8

      # Determine whether we are in AD...
      prod = gee * M0 * T0 * p0 * a0 * Tref * pref * Phiinl * Kinl * pid * pib * pifn * pitn * Gearf
      prod *= pi_fan_des * pi_lpc_des * pi_hpc_des * pi_hpt_des * pi_lpt_des
      prod *= mb_fan_des * mb_lpc_des * mb_hpc_des * mb_hpt_des * mb_lpt_des
      prod *= Nb_fan_des * Nb_lpc_des * Nb_hpc_des * Nb_hpt_des * Nb_lpt_des
      prod *= A2 * A2_5 * A5 * A7
      prod *= Ttf * ifuel * etab
      prod *= epf0 * eplc0 * ephc0 * epht0 * eplt0
      prod *= mofft * Pofft * Tt25off * pt25off * epsl * epsh
      prod *= Mtexit * dTstrk * StA * efilm * tfilm
      prod *= M4a * ruc
      prod *= epsrow[1] * M2 * pif * pilc * pihc * mbf * mblc * mbhc * Tt4 * pt8 * mcore * M2_5
      T = typeof(prod)

      #---- typed turbine component instances (used by turbine_efficiency and turbine_mb_residual)
      turb_hp = Turbine(pi_hpt_des, mb_hpt_des, Nb_hpt_des, epht0;
                        map = TurbineMap{T}(T(Tmaph[1]), T(Tmaph[2])))
      turb_lp = Turbine(pi_lpt_des, mb_lpt_des, Nb_lpt_des, eplt0;
                        map = TurbineMap{T}(T(Tmapl[1]), T(Tmapl[2])))

      #---- Newton system arrays
      res = zeros(T, 9, 1)
      a = zeros(T, 9, 9)
      rrel = zeros(T, 9)
      rsav = zeros(T, 9)
      asav = zeros(T, 9, 10)

      res_dlls = zeros(T, 9)
      a_dlls = zeros(T, 9, 9)

      Random.seed!(1234) #Seed for RNG in relaxation

      #---- number of gas constituents
      n = 6

      #---- mass fractions
      alpha = zeros(T, n)    # air
      beta = zeros(T, n)     # fuel
      gamma = zeros(T, n)    # combustion-caused change in air
      lambda = zeros(T, n)   # combustion product gas
      lambdap = zeros(T, n)  # combustion product gas with cooling flow added

      lam_Tt3 = zeros(T, n)
      lam_Ttf = zeros(T, n)
      lam_pl = zeros(T, n)
      lam_ph = zeros(T, n)
      lam_ml = zeros(T, n)
      lam_mh = zeros(T, n)
      lam_Tb = zeros(T, n)
      lamp_pl = zeros(T, n)
      lamp_ph = zeros(T, n)
      lamp_mf = zeros(T, n)
      lamp_ml = zeros(T, n)
      lamp_mh = zeros(T, n)
      lamp_Tb = zeros(T, n)
      lamp_Mi = zeros(T, n)

      T_al = zeros(T, n)
      p_al = zeros(T, n)
      h_al = zeros(T, n)
      s_al = zeros(T, n)
      R_al = zeros(T, n)
      cp_al = zeros(T, n)

      #---- minimum allowable fan efficiency
      epfmin = 0.60

      #---- typed compressor component instances (used by compressor_pratd, which calls compressor_efficiency internally)
      comp_fan = Compressor(pi_fan_des,  mb_fan_des,  Nb_fan_des,  epf0,  T(0.60), FanMap, windmilling=true)
      comp_lpc = Compressor(pi_lpc_des, mb_lpc_des, Nb_lpc_des, eplc0, T(0.70), LPCMap)
      comp_hpc = Compressor(pi_hpc_des, mb_hpc_des, Nb_hpc_des, ephc0, T(0.70), HPCMap)

      #---- typed shaft component instances (used by hp_shaft_workd, lp_shaft_workd, shaft_speed_residual)
      shaft_hp = Shaft(epsh, T(1.0))
      shaft_lp = Shaft(epsl, Gearf)

      #---- splitter component instance (used by bypass_ratio)
      splitter = Splitter()

      #---- typed combustor component instance (used by combustor_burnd)
      burner = Combustor(pib, etab, Ttf, ifuel, hvap)

      #---- nozzle component instances (used by nozzle_exit and nozzle_massflow_residual)
      #     pn = 1 since pt18 = pt13*pifn and pt8 = pt5c*pitn already incorporate the loss
      nozzle_fan  = Nozzle(one(T), T(A7))
      nozzle_core = Nozzle(one(T), T(A5))

      #---- inlet component instance (pid used in diffuser; Kinl/Phiinl/eng_has_BLI_cores for BLI)
      inlet = Inlet(pid; Kinl=Kinl, Phiinl=Phiinl, eng_has_BLI_cores=eng_has_BLI_cores)

      # from "airfrac.inc"
      # air fractions  
      #        N2      O2      CO2    H2O      Ar       fuel
      alpha = [AIR_ALPHA..., 0.0]
      beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

      #---- convergence tolerance
      #      data toler  1.0e-7 
      toler = 1.0e-9

      itmax = 50

      #---- max fan-face Mach number, above which it will be artificially limited
      Mimax = 0.98

      Lprint = false

      if opt_calc_call == CalcMode.FixedTt4OffDes
            Tt4spec = Tt4
      elseif opt_calc_call == CalcMode.FixedFeOffDes
            Fspec = Feng
      end

      #---- number of air constitutents (all but fuel)
      nair = n - 1

      # ===============================================================
      #---- set combustion-change mass fractions gamma[i] for specified fuel
      gamma = gasfuel(ifuel, n)

      # Convert type for ForwardDiff
      if (typeof(etab) <: ForwardDiff.Dual)
            gamma = convert(Array{typeof(etab)}, gamma)
      end
      #---- apply combustor efficiency
      for i = 1:nair
            gamma[i] = etab * gamma[i]
      end
      gamma[n] = 1.0 - etab

      #
      # ===============================================================
      #---- freestream static quantities
      s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0)
      gam0 = cp0 / (cp0 - R0)
      #      a0 = sqrt(gam0*R0*T0)
      u0 = M0 * a0
      #
      # ===============================================================
      #---- freestream total quantities
      hspec = h0 + 0.5 * u0^2
      Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
      Tt0 = gas_tset(alpha, nair, hspec, Tguess)
      st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
      pt0 = p0 * exp((st0 - s0) / Rt0)
      at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))

      #
      # ===============================================================
      #---- Offtake plume flow 9
      Trat = (p0 / pt25off)^(Rt0 / cpt0)
      if (Trat < 1.0)
            u9 = sqrt(2.0 * cpt0 * Tt25off * (1.0 - Trat))
            rho9 = p0 / (Rt0 * Tt0 * Trat)
      else
            u9 = 0.0
            rho9 = p0 / (Rt0 * Tt0)
      end

      # ===============================================================
      #---- diffuser flow 0-2
      Tt12 = Tt0
      st12 = st0
      ht12 = ht0
      cpt12 = cpt0
      Rt12 = Rt0
      pt12 = pt0 * inlet.pid

      #---- initial guesses for primary Newton variables
      pf = pif
      pl = pilc
      ph = pihc
      mf = mbf
      ml = mblc
      mh = mbhc
      Tb = Tt4
      Pc = pt8
      Mi = M2

      if (pf == 0.0)
            pf = pi_fan_des
      end
      if (pl == 0.0)
            pl = pi_lpc_des
      end
      if (ph == 0.0)
            ph = pi_hpc_des
      end
      if (mf == 0.0)
            mf = mb_fan_des
      end
      if (ml == 0.0)
            ml = mb_lpc_des
      end
      if (mh == 0.0)
            mh = mb_hpc_des
      end
      if (Mi == 0.0)
            Mi = 0.6
      end
      if (Mi == 1.0)
            Mi = 0.6
      end

      #---- total cooling mass flow ratio
      fc = 0.0

      if opt_cooling == CoolingOpt.FixedCoolingFlowRatio
            if (mcore == 0.0)
                  fo = 0.0
            else
                  fo = mofft / mcore
            end
            for icrow = 1:ncrow
                  fc = fc + (1.0 - fo) * epsrow[icrow]
            end
      end

      fc_pl = 0.0
      fc_ph = 0.0
      fc_ml = 0.0
      fc_mh = 0.0
      fc_Tb = 0.0

      for iter = 1:itmax
      
            if (iter == -1)
                  eps1 = 2.0e-7
                  eps = eps1

                  pf = pf + eps
                  j = 1

            end

            epf = epf0
            eplc = eplc0
            ephc = ephc0
            epht = epht0
            eplt = eplt0

            # println("epf, eplc, ephc, epht, eplt")
            # println(epf, eplc, ephc, epht, eplt)
            # exit()

            # HSC: SEEMS TO BE FINE

            # ===============================================================
            #---- set fan inlet total pressure pt13 corrected for BLI
            #-    (Tt2,pt2 approximated with Tt0,pt0 here to avoid circular definition)
            if (u0 == 0.0)
                  #----- static case... no BL defect
                  sbfan = 0.0
                  sbfan_mf = 0.0
                  sbfan_ml = 0.0
                  sbfan_Mi = 0.0

                  sbcore = 0.0
                  sbcore_mf = 0.0
                  sbcore_ml = 0.0
                  sbcore_Mi = 0.0

            else
                  #----- account for inlet BLI defect via mass-averaged entropy
                  a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2)
                  a2sq_Mi = -a2sq / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2) *
                            (gam0 - 1.0) * Mi

                  if eng_has_BLI_cores
                        #------ BL mixes with fan + core flow
                        #c      mmix    = mf*sqrt(Tref/Tt2 ) * pt2    /pref
                        #c             + ml*sqrt(Tref/Tt2a) * pt2a   /pref
                        #c      mmix_mf =    sqrt(Tref/Tt2 ) * pt2    /pref
                        #c      mmix_ml =    sqrt(Tref/Tt2a) * pt2a   /pref
                        #c      mmix_Mi = mf*sqrt(Tref/Tt2 ) * pt2_Mi /pref
                        #c             + ml*sqrt(Tref/Tt2a) * pt1_9_Mi/pref

                        mmix = mf * sqrt(Tref / Tt0) * pt0 / pref +
                               ml * sqrt(Tref / Tt0) * pt0 / pref
                        mmix_mf = sqrt(Tref / Tt0) * pt0 / pref
                        mmix_ml = sqrt(Tref / Tt0) * pt0 / pref
                        mmix_Mi = 0.0

                        sbfan = Kinl * gam0 / (mmix * a2sq)
                        sbfan_mf = (-sbfan / mmix) * mmix_mf
                        sbfan_ml = (-sbfan / mmix) * mmix_ml
                        sbfan_Mi = (-sbfan / mmix) * mmix_Mi +
                                   (-sbfan / a2sq) * a2sq_Mi

                        sbcore = sbfan
                        sbcore_mf = sbfan_mf
                        sbcore_ml = sbfan_ml
                        sbcore_Mi = sbfan_Mi

                  else #clean flow
                        #------ BL mixes with fan flow only
                        #c      mmix    = mf*sqrt(Tref/Tt2) * pt2   /pref
                        #c      mmix_mf =    sqrt(Tref/Tt2) * pt2   /pref
                        #c      mmix_Mi = mf*sqrt(Tref/Tt2) * pt2_Mi/pref

                        mmix = mf * sqrt(Tref / Tt0) * pt0 / pref
                        mmix_mf = sqrt(Tref / Tt0) * pt0 / pref
                        mmix_Mi = 0.0

                        sbfan = Kinl * gam0 / (mmix * a2sq)
                        sbfan_mf = (-sbfan / mmix) * mmix_mf
                        sbfan_ml = 0.0
                        sbfan_Mi = (-sbfan / mmix) * mmix_Mi +
                                    (-sbfan / a2sq) * a2sq_Mi

                        sbcore = 0.0
                        sbcore_mf = 0.0
                        sbcore_ml = 0.0
                        sbcore_Mi = 0.0
                  end
            end

            #---- note: BL is assumed adiabatic, 
            #-     so Tt2,ht2,st2,cpt2,Rt2  will not change due to BL ingestion
            Tt2 = Tt12
            ht2 = ht12
            st2 = st12
            cpt2 = cpt12
            Rt2 = Rt12
            pt2 = pt12 * exp(-sbfan)
            pt2_mf = pt2 * (-sbfan_mf)
            pt2_ml = pt2 * (-sbfan_ml)
            pt2_Mi = pt2 * (-sbfan_Mi)

            Tt2a = Tt12
            ht2a = ht12
            st2a = st12
            cpt2a = cpt12
            Tt1_9_ht1_9 = 1 / cpt2a
            Rt2a = Rt12
            pt2a = pt12 * exp(-sbcore)
            pt1_9_mf = pt2a * (-sbcore_mf)
            pt1_9_ml = pt2a * (-sbcore_ml)
            pt1_9_Mi = pt2a * (-sbcore_Mi)

            p2, T2, h2, s2, cp2, R2,
            p2_st2,
            p2_pt2,
            dum,
            p2_Tt2, T2_Tt2, h2_Tt2, s2_Tt2,
            p2_ht2, T2_ht2, h2_ht2, s2_ht2,
            p2_Mi, T2_Mi, h2_Mi, s2_Mi,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_machd(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, Mi, 1.0)
            u2 = sqrt(2.0 * (ht2 - h2))
            u2_Mi = (-1.0 / u2) * h2_Mi

            rho2 = p2 / (R2 * T2)
            rho2_p2 = 1.0 / (R2 * T2)
            rho2_T2 = -rho2 / T2

            rho2_Mi = rho2_T2 * T2_Mi +
                      rho2_p2 * (p2_pt2 * pt2_Mi + p2_Mi)
            rho2_mf = rho2_p2 * p2_pt2 * pt2_mf
            rho2_ml = rho2_p2 * p2_pt2 * pt2_ml


            p1_9, T1_9, h1_9, s1_9, cp1_9, R1_9,
            p1_9_st1_9,
            p1_9_pt1_9,
            dum,
            p1_9_Tt1_9, T1_9_Tt1_9, h1_9_Tt1_9, s1_9_Tt1_9,
            p1_9_ht1_9, T1_9_ht1_9, h1_9_ht1_9, s1_9_ht1_9,
            p1_9_Mi, T1_9_Mi, h1_9_Mi, s1_9_Mi,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_machd(alpha, nair,
                  pt2a, Tt2a, ht2a, st2a, cpt2a, Rt2a, 0.0, Mi, 1.0)
            u1_9 = sqrt(2.0 * (ht2a - h1_9))
            u1_9_Mi = (-1.0 / u1_9) * h1_9_Mi

            rho1_9 = p1_9 / (R1_9 * T1_9)
            rho1_9_p1_9 = 1.0 / (R1_9 * T1_9)
            rho1_9_T1_9 = -rho1_9 / T1_9

            rho1_9_Mi = rho1_9_T1_9 * T1_9_Mi +
                       rho1_9_p1_9 * (p1_9_pt1_9 * pt1_9_Mi + p1_9_Mi)
            rho1_9_mf = rho1_9_p1_9 * p1_9_pt1_9 * pt1_9_mf
            rho1_9_ml = rho1_9_p1_9 * p1_9_pt1_9 * pt1_9_ml

            # ===============================================================
            #---- fan flow 2-7
            pt13, Tt13, ht13, st13, cpt13, Rt13,
            pt2_1_pt2,
            pt2_1_st2, Tt2_1_st2, ht2_1_st2, st2_1_st2,
            pt2_1_pf, Tt2_1_pf, ht2_1_pf, st2_1_pf,
            pt2_1_mf_eff, Tt2_1_mf_eff, ht2_1_mf_eff, st2_1_mf_eff,
            Nf, Nf_pf, Nf_mf, epf, epf_pf, epf_mf =
                  compressor_pratd(comp_fan, alpha, nair,
                        pt2, Tt2, ht2, st2, cpt2, Rt2, pf, mf)

            pt2_1_mf = pt2_1_mf_eff + pt2_1_pt2 * pt2_mf
            Tt2_1_mf = Tt2_1_mf_eff
            ht2_1_mf = ht2_1_mf_eff
            st2_1_mf = st2_1_mf_eff

            pt2_1_ml = pt2_1_pt2 * pt2_ml
            pt2_1_Mi = pt2_1_pt2 * pt2_Mi

            #---- fan duct nozzle total quantities
            pt18 = pt13 * pifn
            Tt18 = Tt13
            ht18 = ht13
            st18 = st13
            cpt18 = cpt13
            Rt18 = Rt13

            pt7_pf = pt2_1_pf * pifn
            Tt7_pf = Tt2_1_pf
            ht7_pf = ht2_1_pf
            st7_pf = st2_1_pf

            pt7_mf = pt2_1_mf * pifn
            Tt7_mf = Tt2_1_mf
            ht7_mf = ht2_1_mf
            st7_mf = st2_1_mf

            pt7_ml = pt2_1_ml * pifn
            ht7_ml = 0.0

            pt7_Mi = pt2_1_Mi * pifn
            ht7_Mi = 0.0
            # ===============================================================
            #---- Compressor precooler 19-19c
            pt2ac = pt2a - Δp_PreC
            ht2ac = ht2a + Δh_PreC
            Tt2ac, Tt1_9c_ht1_9c, _ = gas_tsetd(alpha, nair, ht2ac, Tt2a)
            st2ac, st1_9c_Tt1_9c, ht2ac, ht29c_Tt29c, cpt2ac, cpt1_9c_Tt1_9c, Rt2ac = gassumd(alpha, nair, Tt2ac)

             #Derivatives with respect to pressure are unchanged because pt1_9c_pt1_9 = 1
            pt1_9c_pt1_9 = 1.0
            ht1_9c_ht1_9 = 1.0

            Tt1_9c_Tt1_9 = Tt1_9c_ht1_9c * ht1_9c_ht1_9 / Tt1_9_ht1_9
            pt1_9c_mf = pt1_9_mf * pt1_9c_pt1_9
            pt1_9c_Mi = pt1_9_Mi * pt1_9c_pt1_9
            pt1_9c_ml = pt1_9_ml * pt1_9c_pt1_9

            p1_9c, T1_9c, h1_9c, s1_9c, cp1_9c, R1_9c,
            p1_9c_st1_9c,
            p1_9c_pt1_9c,
            dum,
            p1_9c_Tt1_9c, T1_9c_Tt1_9c, h1_9c_Tt1_9c, s1_9c_Tt1_9c,
            p1_9c_ht1_9c, T1_9c_ht1_9c, h1_9c_ht1_9c, s1_9c_ht1_9c,
            p1_9c_Mi, T1_9c_Mi, h1_9c_Mi, s1_9c_Mi,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_machd(alpha, nair,
                  pt2ac, Tt2ac, ht2ac, st2ac, cpt2ac, Rt2ac, 0.0, Mi, 1.0)
            u1_9c = sqrt(2.0 * (ht2ac - h1_9c))
            u1_9c_Mi = (-1.0 / u1_9c) * h1_9c_Mi

            rho1_9c = p1_9c / (R1_9c * T1_9c)
            rho1_9c_p1_9c = 1.0 / (R1_9c * T1_9c)
            rho1_9c_T1_9c = -rho1_9c / T1_9c

            rho1_9c_Mi = rho1_9c_T1_9c * T1_9c_Mi +
                       rho1_9c_p1_9c * (p1_9c_pt1_9c * pt1_9c_Mi + p1_9c_Mi)
            rho1_9c_mf = rho1_9c_p1_9c * p1_9c_pt1_9c * pt1_9c_mf
            rho1_9c_ml = rho1_9c_p1_9c * p1_9c_pt1_9c * pt1_9c_ml


            #--------------------------------------------------------------
            #---- offtake mass ratio
            #cc   fo = mofft / mcore
            fo = mofft / ml * sqrt(Tt2ac / Tref) * pref / pt2ac
            fo_ml = -fo / ml - (fo / pt2ac) * pt1_9c_ml
            fo_mf = -(fo / pt2ac) * pt1_9c_mf
            fo_Mi = -(fo / pt2ac) * pt1_9c_Mi

            #---- normalized power offtake Pofft / mcore
            Pom = Pofft / ml * sqrt(Tt2ac / Tref) * pref / pt2ac
            Pom_ml = -Pom / ml - (Pom / pt2ac) * pt1_9c_ml
            Pom_mf = -(Pom / pt2ac) * pt1_9c_mf
            Pom_Mi = -(Pom / pt2ac) * pt1_9c_Mi
      
            # ===============================================================
            #---- LP compressor flow 2-25
            pt2_5, Tt2_5, ht2_5, st2_5, cpt2_5, Rt2_5,
            pt2_5_pt1_9c,
            pt2_5_st1_9c, Tt2_5_st1_9c, ht2_5_st1_9c, st2_5_st1_9c,
            pt2_5_pl, Tt2_5_pl, ht2_5_pl, st2_5_pl,
            pt2_5_ml_eff, Tt2_5_ml_eff, ht2_5_ml_eff, st2_5_ml_eff,
            Nl, Nl_pl, Nl_ml, eplc, eplc_pl, eplc_ml =
                  compressor_pratd(comp_lpc, alpha, nair,
                        pt2ac, Tt2ac, ht2ac, st2ac, cpt2ac, Rt2ac, pl, ml)

            pt2_5_ml = pt2_5_ml_eff + pt2_5_pt1_9c * pt1_9c_ml
            Tt2_5_ml = Tt2_5_ml_eff
            ht2_5_ml = ht2_5_ml_eff
            st2_5_ml = st2_5_ml_eff

            pt2_5_mf = pt2_5_pt1_9c * pt1_9c_mf
            pt2_5_Mi = pt2_5_pt1_9c * pt1_9c_Mi

            Tt2_5_ht2_5 = 1 / cpt2_5
            st2_5_Tt2_5 = cpt2_5 / Tt2_5

            # ===============================================================
            #---- Compressor intercooler 25-25c
            pt2_5c = pt2_5 - Δp_InterC
            ht2_5c = ht2_5 + Δh_InterC
            Tt2_5c, Tt2_5c_ht2_5c, _ = gas_tsetd(alpha, nair, ht2_5c, Tt2_5)
            st2_5c, st2_5c_Tt2_5c, ht2_5c, ht2_5c_Tt2_5c, cpt2_5c, cpt2_5c_Tt2_5c, Rt2_5c = gassumd(alpha, nair, Tt2_5c)

            #Derivatives with respect to pressure are unchanged
            pt2_5c_pt2_5 = 1.0
            ht2_5c_ht2_5 = 1.0

            Tt2_5c_Tt2_5 = Tt2_5c_ht2_5c * ht2_5c_ht2_5 / Tt2_5_ht2_5
            st2_5c_st2_5 = st2_5c_Tt2_5c * Tt2_5c_Tt2_5 / st2_5_Tt2_5

            pt2_5c_pl = pt2_5_pl * pt2_5c_pt2_5
            Tt2_5c_pl = Tt2_5_pl * Tt2_5c_Tt2_5
            ht2_5c_pl = ht2_5_pl * ht2_5c_ht2_5
            st2_5c_pl = st2_5_pl * st2_5c_st2_5

            pt2_5c_ml = pt2_5_ml * pt2_5c_pt2_5
            Tt2_5c_ml = Tt2_5_ml * Tt2_5c_Tt2_5
            ht2_5c_ml = ht2_5_ml * ht2_5c_ht2_5
            st2_5c_ml = st2_5_ml * st2_5c_st2_5

            pt2_5c_mf = pt2_5_mf * pt2_5c_pt2_5
            pt2_5c_Mi = pt2_5_Mi * pt2_5c_pt2_5
            pt2_5c_ml = pt2_5_ml * pt2_5c_pt2_5
      
            # ===============================================================
            #---- HP compressor flow 25-3
            pt3, Tt3, ht3, st3, cpt3, Rt3,
            pt3_pt2_5c,
            pt3_st2_5c, Tt3_st2_5c, ht3_st2_5c, st3_st2_5c,
            pt3_ph, Tt3_ph, ht3_ph, st3_ph,
            pt3_mh_eff, Tt3_mh_eff, ht3_mh_eff, st3_mh_eff,
            Nh, Nh_ph, Nh_mh, ephc, ephc_ph, ephc_mh =
                  compressor_pratd(comp_hpc, alpha, nair,
                        pt2_5c, Tt2_5c, ht2_5c, st2_5c, cpt2_5c, Rt2_5c, ph, mh)

            pt3_pl = pt3_pt2_5c * pt2_5c_pl +
                     pt3_st2_5c * st2_5c_pl
            Tt3_pl = Tt3_st2_5c * st2_5c_pl
            ht3_pl = ht3_st2_5c * st2_5c_pl
            st3_pl = st3_st2_5c * st2_5c_pl

            pt3_ml = pt3_pt2_5c * pt2_5c_ml +
                     pt3_st2_5c * st2_5c_ml
            Tt3_ml = Tt3_st2_5c * st2_5c_ml
            ht3_ml = ht3_st2_5c * st2_5c_ml
            st3_ml = st3_st2_5c * st2_5c_ml

            pt3_mf = pt3_pt2_5c * pt2_5c_mf

            pt3_mh = pt3_mh_eff
            Tt3_mh = Tt3_mh_eff
            ht3_mh = ht3_mh_eff
            st3_mh = st3_mh_eff

            pt3_Mi = pt3_pt2_5c * pt2_5c_Mi

            Tt3_ht3 = 1 / cpt3

            # ===============================================================
            #---- burner: fuel fraction, exit composition + station 4 state with Jacobian
            #     (combustor_burnd wraps gas_burnd + gassumd + composition chain)
            ffb, lambda,
            Tt4, ht4, st4, cpt4, Rt4,
            ffb_Tt3, ffb_Tb,
            ht4_Tt3, st4_Tt3, cpt4_Tt3, Rt4_Tt3,
            ht4_Tb, st4_Tb, cpt4_Tb, Rt4_Tb,
            lam_Tt3, lam_Tb = combustor_burnd(burner, alpha, nair, Tt3, Tb)
            Tt4_Tb = 1.0  # Tt4 = Tb, needed downstream (cooling, Newton row 7)

            #---- chain combustor derivatives to Newton variables via Tt3
            ffb_pl = ffb_Tt3 * Tt3_pl
            ffb_ph = ffb_Tt3 * Tt3_ph
            ffb_ml = ffb_Tt3 * Tt3_ml
            ffb_mh = ffb_Tt3 * Tt3_mh
            for i = 1:nair
                  lam_pl[i] = lam_Tt3[i] * Tt3_pl
                  lam_ph[i] = lam_Tt3[i] * Tt3_ph
                  lam_ml[i] = lam_Tt3[i] * Tt3_ml
                  lam_mh[i] = lam_Tt3[i] * Tt3_mh
            end

            #---- station 4 state derivatives via Tt3 chain
            st4_pl  = st4_Tt3  * Tt3_pl;  st4_ph  = st4_Tt3  * Tt3_ph
            st4_ml  = st4_Tt3  * Tt3_ml;  st4_mh  = st4_Tt3  * Tt3_mh
            ht4_pl  = ht4_Tt3  * Tt3_pl;  ht4_ph  = ht4_Tt3  * Tt3_ph
            ht4_ml  = ht4_Tt3  * Tt3_ml;  ht4_mh  = ht4_Tt3  * Tt3_mh
            cpt4_pl = cpt4_Tt3 * Tt3_pl;  cpt4_ph = cpt4_Tt3 * Tt3_ph
            cpt4_ml = cpt4_Tt3 * Tt3_ml;  cpt4_mh = cpt4_Tt3 * Tt3_mh
            Rt4_pl  = Rt4_Tt3  * Tt3_pl;  Rt4_ph  = Rt4_Tt3  * Tt3_ph
            Rt4_ml  = Rt4_Tt3  * Tt3_ml;  Rt4_mh  = Rt4_Tt3  * Tt3_mh

            #---- station 4 pressure (depends on pt3, not temperatures)
            pt4 = pib * pt3
            pt4_pl = pib * pt3_pl
            pt4_ph = pib * pt3_ph
            pt4_ml = pib * pt3_ml
            pt4_mh = pib * pt3_mh
            pt4_Tb = 0.0
            pt4_Mi = pib * pt3_Mi

            # HSC: SEEMS TO BE FINE
            # ===============================================================
            if opt_cooling == CoolingOpt.NoCooling
                  #----- no cooling air present... station 4.1 is same as 4
                  pt4_1 = pt4
                  Tt4_1 = Tt4
                  ht4_1 = ht4
                  st4_1 = st4
                  cpt4_1 = cpt4
                  Rt4_1 = Rt4
                  for i = 1:nair
                        lambdap[i] = lambda[i]
                        lamp_pl[i] = lam_pl[i]
                        lamp_ph[i] = lam_ph[i]
                        lamp_mf[i] = 0.0
                        lamp_ml[i] = lam_ml[i]
                        lamp_mh[i] = lam_mh[i]
                        lamp_Tb[i] = lam_Tb[i]
                        lamp_Mi[i] = 0.0
                  end

                  pt4_1_pl = pt4_pl
                  pt4_1_ph = pt4_ph
                  pt4_1_ml = pt4_ml
                  pt4_1_mh = pt4_mh
                  pt4_1_Tb = pt4_Tb
                  pt4_1_Mi = pt4_Mi

                  Tt4_1_pl = 0.0
                  Tt4_1_ph = 0.0
                  Tt4_1_mf = 0.0
                  Tt4_1_ml = 0.0
                  Tt4_1_mh = 0.0
                  Tt4_1_Tb = 1.0
                  Tt4_1_Mi = 0.0

                  ht4_1_pl = ht4_pl
                  ht4_1_ph = ht4_ph
                  ht4_1_mf = 0.0
                  ht4_1_ml = ht4_ml
                  ht4_1_mh = ht4_mh
                  ht4_1_Tb = ht4_Tb
                  ht4_1_Mi = 0.0

                  st4_1_pl = st4_pl
                  st4_1_ph = st4_ph
                  st4_1_ml = st4_ml
                  st4_1_mh = st4_mh
                  st4_1_Tb = st4_Tb

                  cpt4_1_pl = cpt4_pl
                  cpt4_1_ph = cpt4_ph
                  cpt4_1_ml = cpt4_ml
                  cpt4_1_mh = cpt4_mh
                  cpt4_1_Tb = cpt4_Tb

                  Rt4_1_pl = Rt4_pl
                  Rt4_1_ph = Rt4_ph
                  Rt4_1_ml = Rt4_ml
                  Rt4_1_mh = Rt4_mh
                  Rt4_1_Tb = Rt4_Tb
                  #
                  fc = 0.0
                  fc_pl = 0.0
                  fc_ph = 0.0
                  fc_mf = 0.0
                  fc_ml = 0.0
                  fc_mh = 0.0
                  fc_Tb = 0.0
                  fc_Mi = 0.0

                  #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
                  ff = (1.0 - fo) * ffb
                  ff_pl = (1.0 - fo) * ffb_pl
                  ff_ph = (1.0 - fo) * ffb_ph
                  ff_mf = -fo_mf * ffb
                  ff_ml = (1.0 - fo) * ffb_ml - fo_ml * ffb
                  ff_mh = (1.0 - fo) * ffb_mh
                  ff_Tb = (1.0 - fo) * ffb_Tb
                  ff_Mi = -fo_Mi * ffb

                  #----------------------------------------------------------------
            else
                  #----- cooling air is present... 

                  #----- hot-section temperature ratio for each blade row (for cooling model)
                  gmi4 = Rt4 / (cpt4 - Rt4)
                  gmi4_Rt4 = (1.0 + gmi4) / (cpt4 - Rt4)
                  gmi4_cpt4 = -gmi4 / (cpt4 - Rt4)
                  Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)
                  Trr_gmi4 = -Trrat / (1.0 + 0.5 * gmi4 * Mtexit^2) * 0.5 * Mtexit^2
                  Trr_pl = Trr_gmi4 * (gmi4_Rt4 * Rt4_pl + gmi4_cpt4 * cpt4_pl)
                  Trr_ph = Trr_gmi4 * (gmi4_Rt4 * Rt4_ph + gmi4_cpt4 * cpt4_ph)
                  Trr_ml = Trr_gmi4 * (gmi4_Rt4 * Rt4_ml + gmi4_cpt4 * cpt4_ml)
                  Trr_mh = Trr_gmi4 * (gmi4_Rt4 * Rt4_mh + gmi4_cpt4 * cpt4_mh)
                  Trr_Tb = Trr_gmi4 * (gmi4_Rt4 * Rt4_Tb + gmi4_cpt4 * cpt4_Tb)

                  # Heat exchanger to cool turbine cooling air
                  ht_tc = ht3 + Δh_TurbC #Specific enthalpy of turbine cooling air
                  Tt_tc, Tttc_httc, _ = gas_tsetd(alpha, nair, ht_tc, Tt3) #Temperature of turbine cooling air

                  httc_ht3 = 1.0
                  Tttc_Tt3 = Tttc_httc * httc_ht3 / Tt3_ht3

                  if opt_cooling == CoolingOpt.FixedCoolingFlowRatio
                        #------ epsrow(.) is assumed to be passed in.. calculate Tmrow(.)
                        Tmrow_copy = Tmcalc(ncrowx, ncrow,
                              Tt_tc, Tb, dTstrk, Trrat,
                              efilm, tfilm, StA, epsrow)
                        Tmrow[:] = Tmrow_copy[:]

                        #------ total cooling flow fraction
                        fc = 0.0
                        fc_fo = 0.0
                        for icrow = 1:ncrow
                              fc = fc + (1.0 - fo) * epsrow[icrow]
                              fc_fo = fc_fo - epsrow[icrow]
                        end
                        fc_pl = 0.0
                        fc_ph = 0.0
                        fc_mf = fc_fo * fo_mf
                        fc_ml = fc_fo * fo_ml
                        fc_mh = 0.0
                        fc_Tb = 0.0
                        fc_Mi = fc_fo * fo_Mi

                  else
                        #------ calculate cooling mass flow ratios epsrow(.) to get specified Tmrow(.)
                        ncrow, epsrow_copy, epsrow_Tttc, epsrow_Tb, epsrow_Trr = mcool(ncrowx,
                              Tmrow, Tt_tc, Tb, dTstrk, Trrat,
                              efilm, tfilm, StA)
                        epsrow[:] = epsrow_copy[:]

                        #------ total cooling flow fraction
                        fc = 0.0
                        fc_fo = 0.0
                        fc_Tttc = 0.0
                        fc_Tb = 0.0
                        fc_Trr = 0.0
                        for icrow = 1:ncrow
                              fc = fc + (1.0 - fo) * epsrow[icrow]
                              fc_fo = fc_fo - epsrow[icrow]
                              fc_Tttc = fc_Tttc + (1.0 - fo) * epsrow_Tttc[icrow]
                              fc_Tb = fc_Tb + (1.0 - fo) * epsrow_Tb[icrow]
                              fc_Trr = fc_Trr + (1.0 - fo) * epsrow_Trr[icrow]
                        end
                        fc_Tt3 = fc_Tttc * Tttc_Tt3

                        fc_pl = fc_Tt3 * Tt3_pl + fc_Trr * Trr_pl
                        fc_ph = fc_Tt3 * Tt3_ph + fc_Trr * Trr_ph
                        fc_mf = fc_fo * fo_mf
                        fc_ml = fc_Tt3 * Tt3_ml + fc_Trr * Trr_ml + fc_fo * fo_ml
                        fc_mh = fc_Tt3 * Tt3_mh + fc_Trr * Trr_mh
                        fc_Tb = fc_Tb + fc_Trr * Trr_Tb
                        fc_Mi = fc_fo * fo_Mi

                  end

                  #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
                  ff = (1.0 - fo - fc) * ffb
                  ff_pl = (1.0 - fo - fc) * ffb_pl - fc_pl * ffb
                  ff_ph = (1.0 - fo - fc) * ffb_ph - fc_ph * ffb
                  ff_mf = -fc_mf * ffb - fo_mf * ffb
                  ff_ml = (1.0 - fo - fc) * ffb_ml - fc_ml * ffb - fo_ml * ffb
                  ff_mh = (1.0 - fo - fc) * ffb_mh - fc_mh * ffb
                  ff_Tb = (1.0 - fo - fc) * ffb_Tb - fc_Tb * ffb
                  ff_Mi = -fc_Mi * ffb - fo_Mi * ffb

                  #----- calculate station 4.1
                  pt4a = pt4
                  Tt4a = Tt4
                  ht4a = ht4
                  st4a = st4
                  cpt4a = cpt4
                  rt4a = Rt4

                  Tt4a_Tb = Tt4_Tb

                  pt4a_pl = pt4_pl
                  pt4a_ph = pt4_ph
                  pt4a_ml = pt4_ml
                  pt4a_mh = pt4_mh
                  pt4a_Tb = pt4_Tb
                  pt4a_Mi = pt4_Mi

                  ht4a_pl = ht4_pl
                  ht4a_ph = ht4_ph
                  ht4a_ml = ht4_ml
                  ht4a_mh = ht4_mh
                  ht4a_Tb = ht4_Tb

                  st4a_pl = st4_pl
                  st4a_ph = st4_ph
                  st4a_ml = st4_ml
                  st4a_mh = st4_mh
                  st4a_Tb = st4_Tb

                  #----- speed at start-of-mixing station 4a
                  p4a, t4a, h4a, s4a, cp4a, r4a,
                  p4a_st4a,
                  p4a_pt4a,
                  dum,
                  p4a_Tt4a, T4a_Tt4a, h4a_Tt4a, s4a_Tt4a,
                  p4a_ht4a, T4a_ht4a, h4a_ht4a, s4a_ht4a,
                  p4a_M4a, T4a_M4a, h4a_M4a, s4a_M4a,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_machd(lambda, nair,
                        pt4a, Tt4a, ht4a, st4a, cpt4a, rt4a, 0.0, M4a, 1.0)

                  p4a_pl = p4a_st4a * st4a_pl +
                           p4a_pt4a * pt4a_pl +
                           p4a_ht4a * ht4a_pl
                  T4a_pl = T4a_ht4a * ht4a_pl
                  h4a_pl = h4a_ht4a * ht4a_pl
                  s4a_pl = s4a_ht4a * ht4a_pl

                  p4a_ph = p4a_st4a * st4a_ph +
                           p4a_pt4a * pt4a_ph +
                           p4a_ht4a * ht4a_ph
                  T4a_ph = T4a_ht4a * ht4a_ph
                  h4a_ph = h4a_ht4a * ht4a_ph
                  s4a_ph = s4a_ht4a * ht4a_ph


                  p4a_ml = p4a_st4a * st4a_ml +
                           p4a_pt4a * pt4a_ml +
                           p4a_ht4a * ht4a_ml
                  T4a_ml = T4a_ht4a * ht4a_ml
                  h4a_ml = h4a_ht4a * ht4a_ml
                  s4a_ml = s4a_ht4a * ht4a_ml

                  p4a_mh = p4a_st4a * st4a_mh +
                           p4a_pt4a * pt4a_mh +
                           p4a_ht4a * ht4a_mh
                  T4a_mh = T4a_ht4a * ht4a_mh
                  h4a_mh = h4a_ht4a * ht4a_mh
                  s4a_mh = s4a_ht4a * ht4a_mh

                  p4a_Tb = p4a_st4a * st4a_Tb +
                           p4a_pt4a * pt4a_Tb +
                           p4a_ht4a * ht4a_Tb + p4a_Tt4a * Tt4a_Tb
                  T4a_Tb = T4a_ht4a * ht4a_Tb + T4a_Tt4a * Tt4a_Tb
                  h4a_Tb = h4a_ht4a * ht4a_Tb + h4a_Tt4a * Tt4a_Tb
                  s4a_Tb = s4a_ht4a * ht4a_Tb + s4a_Tt4a * Tt4a_Tb

                  p4a_Mi = p4a_pt4a * pt4a_Mi

                  for i = 1:nair
                        p4a_pl = p4a_pl + p_al[i] * lam_pl[i]
                        p4a_ph = p4a_ph + p_al[i] * lam_ph[i]
                        p4a_ml = p4a_ml + p_al[i] * lam_ml[i]
                        p4a_mh = p4a_mh + p_al[i] * lam_mh[i]
                        p4a_Tb = p4a_Tb + p_al[i] * lam_Tb[i]

                        T4a_pl = T4a_pl + T_al[i] * lam_pl[i]
                        T4a_ph = T4a_ph + T_al[i] * lam_ph[i]
                        T4a_ml = T4a_ml + T_al[i] * lam_ml[i]
                        T4a_mh = T4a_mh + T_al[i] * lam_mh[i]
                        T4a_Tb = T4a_Tb + T_al[i] * lam_Tb[i]

                        h4a_pl = h4a_pl + h_al[i] * lam_pl[i]
                        h4a_ph = h4a_ph + h_al[i] * lam_ph[i]
                        h4a_ml = h4a_ml + h_al[i] * lam_ml[i]
                        h4a_mh = h4a_mh + h_al[i] * lam_mh[i]
                        h4a_Tb = h4a_Tb + h_al[i] * lam_Tb[i]

                        s4a_pl = s4a_pl + s_al[i] * lam_pl[i]
                        s4a_ph = s4a_ph + s_al[i] * lam_ph[i]
                        s4a_ml = s4a_ml + s_al[i] * lam_ml[i]
                        s4a_mh = s4a_mh + s_al[i] * lam_mh[i]
                        s4a_Tb = s4a_Tb + s_al[i] * lam_Tb[i]
                  end

                  if (ht4a > h4a)
                        u4a = sqrt(2.0 * (ht4a - h4a))
                        u4a_pl = (ht4a_pl - h4a_pl) / u4a
                        u4a_ph = (ht4a_ph - h4a_ph) / u4a
                        u4a_ml = (ht4a_ml - h4a_ml) / u4a
                        u4a_mh = (ht4a_mh - h4a_mh) / u4a
                        u4a_Tb = (ht4a_Tb - h4a_Tb) / u4a
                  else
                        u4a = 0.0
                        u4a_pl = 0.0
                        u4a_ph = 0.0
                        u4a_ml = 0.0
                        u4a_mh = 0.0
                        u4a_Tb = 0.0
                  end

                  #----- exit speed of cooling air at station 4a
                  uc = ruc * u4a
                  uc_pl = ruc * u4a_pl
                  uc_ph = ruc * u4a_ph
                  uc_ml = ruc * u4a_ml
                  uc_mh = ruc * u4a_mh
                  uc_Tb = ruc * u4a_Tb

                  #----- IGV exit mixing
                  frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
                  frac4_fo = -(1.0 - frac4) / (1.0 - fo + ff)
                  frac4_ff = (1.0 - frac4) / (1.0 - fo + ff)
                  frac4_fc = -1.0 / (1.0 - fo + ff)

                  fracm = fc / (1.0 - fo + ff)
                  fracm_fo = fracm / (1.0 - fo + ff)
                  fracm_ff = -fracm / (1.0 - fo + ff)
                  fracm_fc = 1.0 / (1.0 - fo + ff)

                  frac4_pl = frac4_fc * fc_pl + frac4_ff * ff_pl
                  frac4_ph = frac4_fc * fc_ph + frac4_ff * ff_ph
                  frac4_mf = frac4_fc * fc_mf + frac4_ff * ff_mf + frac4_fo * fo_mf
                  frac4_ml = frac4_fc * fc_ml + frac4_ff * ff_ml + frac4_fo * fo_ml
                  frac4_mh = frac4_fc * fc_mh + frac4_ff * ff_mh
                  frac4_Tb = frac4_fc * fc_Tb + frac4_ff * ff_Tb
                  frac4_Mi = frac4_fc * fc_Mi + frac4_ff * ff_Mi + frac4_fo * fo_Mi

                  fracm_pl = fracm_fc * fc_pl + fracm_ff * ff_pl
                  fracm_ph = fracm_fc * fc_ph + fracm_ff * ff_ph
                  fracm_mf = fracm_fc * fc_mf + fracm_ff * ff_mf + fracm_fo * fo_mf
                  fracm_ml = fracm_fc * fc_ml + fracm_ff * ff_ml + fracm_fo * fo_ml
                  fracm_mh = fracm_fc * fc_mh + fracm_ff * ff_mh
                  fracm_Tb = fracm_fc * fc_Tb + fracm_ff * ff_Tb
                  fracm_Mi = fracm_fc * fc_Mi + fracm_ff * ff_Mi + fracm_fo * fo_Mi

                  #----- mixed constituent fraction vector from mass equation
                  for i = 1:nair
                        lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
                        lamp_pl[i] = frac4_pl * lambda[i] + fracm_pl * alpha[i] +
                                     frac4 * lam_pl[i]
                        lamp_ph[i] = frac4_ph * lambda[i] + fracm_ph * alpha[i] +
                                     frac4 * lam_ph[i]
                        lamp_mf[i] = frac4_mf * lambda[i] + fracm_mf * alpha[i]
                        lamp_ml[i] = frac4_ml * lambda[i] + fracm_ml * alpha[i] +
                                     frac4 * lam_ml[i]
                        lamp_mh[i] = frac4_mh * lambda[i] + fracm_mh * alpha[i] +
                                     frac4 * lam_mh[i]
                        lamp_Tb[i] = frac4_Tb * lambda[i] + fracm_Tb * alpha[i] +
                                     frac4 * lam_Tb[i]
                        lamp_Mi[i] = frac4_Mi * lambda[i] + fracm_Mi * alpha[i]
                  end

                  #derivatives for turbine cooling air when there is a HX
                  httc_pl = ht3_pl * httc_ht3
                  httc_ph = ht3_ph * httc_ht3
                  httc_ml = ht3_ml * httc_ht3
                  httc_mh = ht3_mh * httc_ht3

                  #----- mixed total enthalpy from enthalpy equation
                  ht4_1 = frac4 * ht4 + fracm * ht_tc
                  ht4_1_pl = frac4_pl * ht4 + frac4 * ht4_pl + fracm_pl * ht_tc + fracm * httc_pl
                  ht4_1_ph = frac4_ph * ht4 + frac4 * ht4_ph + fracm_ph * ht_tc + fracm * httc_ph
                  ht4_1_mf = frac4_mf * ht4 + fracm_mf * ht_tc
                  ht4_1_ml = frac4_ml * ht4 + frac4 * ht4_ml + fracm_ml * ht_tc + fracm * httc_ml
                  ht4_1_mh = frac4_mh * ht4 + frac4 * ht4_mh + fracm_mh * ht_tc + fracm * httc_mh
                  ht4_1_Tb = frac4_Tb * ht4 + frac4 * ht4_Tb + fracm_Tb * ht_tc
                  ht4_1_Mi = frac4_mf * ht4 + fracm_Mi * ht_tc

                  #----- total temperature from total enthalpy
                  Tguess = frac4 * Tt4 + fracm * Tt_tc
                  Tt4_1, Tt4_1_ht4_1, T_al = gas_tsetd(lambdap, nair, ht4_1, Tguess)
                  Tt4_1_pl = Tt4_1_ht4_1 * ht4_1_pl
                  Tt4_1_ph = Tt4_1_ht4_1 * ht4_1_ph
                  Tt4_1_mf = Tt4_1_ht4_1 * ht4_1_mf
                  Tt4_1_ml = Tt4_1_ht4_1 * ht4_1_ml
                  Tt4_1_mh = Tt4_1_ht4_1 * ht4_1_mh
                  Tt4_1_Tb = Tt4_1_ht4_1 * ht4_1_Tb
                  Tt4_1_Mi = Tt4_1_ht4_1 * ht4_1_Mi
                  for i = 1:nair
                        Tt4_1_pl = Tt4_1_pl + T_al[i] * lamp_pl[i]
                        Tt4_1_ph = Tt4_1_ph + T_al[i] * lamp_ph[i]
                        Tt4_1_mf = Tt4_1_mf + T_al[i] * lamp_mf[i]
                        Tt4_1_ml = Tt4_1_ml + T_al[i] * lamp_ml[i]
                        Tt4_1_mh = Tt4_1_mh + T_al[i] * lamp_mh[i]
                        Tt4_1_Tb = Tt4_1_Tb + T_al[i] * lamp_Tb[i]
                        Tt4_1_Mi = Tt4_1_Mi + T_al[i] * lamp_Mi[i]
                  end

                  #----- will also need st4_1,cpt4_1,Rt4_1 derivatives
                  st4_1, st4_1_Tt4_1,
                  ht4_1, ht4_1_Tt4_1,
                  cpt4_1, cpt4_1_Tt4_1, Rt4_1 = gassumd(lambdap, nair, Tt4_1)
                  st4_1_pl = st4_1_Tt4_1 * Tt4_1_pl
                  st4_1_ph = st4_1_Tt4_1 * Tt4_1_ph
                  st4_1_mf = st4_1_Tt4_1 * Tt4_1_mf
                  st4_1_ml = st4_1_Tt4_1 * Tt4_1_ml
                  st4_1_mh = st4_1_Tt4_1 * Tt4_1_mh
                  st4_1_Tb = st4_1_Tt4_1 * Tt4_1_Tb
                  st4_1_Mi = st4_1_Tt4_1 * Tt4_1_Mi

                  cpt4_1_pl = cpt4_1_Tt4_1 * Tt4_1_pl
                  cpt4_1_ph = cpt4_1_Tt4_1 * Tt4_1_ph
                  cpt4_1_mf = cpt4_1_Tt4_1 * Tt4_1_mf
                  cpt4_1_ml = cpt4_1_Tt4_1 * Tt4_1_ml
                  cpt4_1_mh = cpt4_1_Tt4_1 * Tt4_1_mh
                  cpt4_1_Tb = cpt4_1_Tt4_1 * Tt4_1_Tb
                  cpt4_1_Mi = cpt4_1_Tt4_1 * Tt4_1_Mi

                  Rt4_1_pl = 0.0
                  Rt4_1_ph = 0.0
                  Rt4_1_mf = 0.0
                  Rt4_1_ml = 0.0
                  Rt4_1_mh = 0.0
                  Rt4_1_Tb = 0.0
                  Rt4_1_Mi = 0.0
                  for i = 1:nair
                        si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, Tt4_1)
                        st4_1_pl = st4_1_pl + si * lamp_pl[i]
                        st4_1_ph = st4_1_ph + si * lamp_ph[i]
                        st4_1_mf = st4_1_mf + si * lamp_mf[i]
                        st4_1_ml = st4_1_ml + si * lamp_ml[i]
                        st4_1_mh = st4_1_mh + si * lamp_mh[i]
                        st4_1_Tb = st4_1_Tb + si * lamp_Tb[i]
                        st4_1_Mi = st4_1_Mi + si * lamp_Mi[i]

                        cpt4_1_pl = cpt4_1_pl + cpi * lamp_pl[i]
                        cpt4_1_ph = cpt4_1_ph + cpi * lamp_ph[i]
                        cpt4_1_mf = cpt4_1_mf + cpi * lamp_mf[i]
                        cpt4_1_ml = cpt4_1_ml + cpi * lamp_ml[i]
                        cpt4_1_mh = cpt4_1_mh + cpi * lamp_mh[i]
                        cpt4_1_Tb = cpt4_1_Tb + cpi * lamp_Tb[i]
                        cpt4_1_Mi = cpt4_1_Mi + cpi * lamp_Mi[i]

                        Rt4_1_pl = Rt4_1_pl + Ri * lamp_pl[i]
                        Rt4_1_ph = Rt4_1_ph + Ri * lamp_ph[i]
                        Rt4_1_mf = Rt4_1_mf + Ri * lamp_mf[i]
                        Rt4_1_ml = Rt4_1_ml + Ri * lamp_ml[i]
                        Rt4_1_mh = Rt4_1_mh + Ri * lamp_mh[i]
                        Rt4_1_Tb = Rt4_1_Tb + Ri * lamp_Tb[i]
                        Rt4_1_Mi = Rt4_1_Mi + Ri * lamp_Mi[i]
                  end

                  #----- mixed velocity from momentum equation, assuming constant static pressure
                  p4_1 = p4a
                  p4_1_pl = p4a_pl
                  p4_1_ph = p4a_ph
                  p4_1_ml = p4a_ml
                  p4_1_mh = p4a_mh
                  p4_1_Tb = p4a_Tb
                  p4_1_Mi = p4a_Mi

                  u4_1 = frac4 * u4a + fracm * uc
                  u4_1_pl = frac4_pl * u4a + frac4 * u4a_pl + fracm_pl * uc + fracm * uc_pl
                  u4_1_ph = frac4_ph * u4a + frac4 * u4a_ph + fracm_ph * uc + fracm * uc_ph
                  u4_1_mf = frac4_mf * u4a + fracm_mf * uc
                  u4_1_ml = frac4_ml * u4a + frac4 * u4a_ml + fracm_ml * uc + fracm * uc_ml
                  u4_1_mh = frac4_mh * u4a + frac4 * u4a_mh + fracm_mh * uc + fracm * uc_mh
                  u4_1_Tb = frac4_Tb * u4a + frac4 * u4a_Tb + fracm_Tb * uc + fracm * uc_Tb
                  u4_1_Mi = frac4_Mi * u4a + fracm_Mi * uc

                  #----- static temperature from static enthalpy
                  h4_1 = ht4_1 - 0.5 * u4_1^2
                  h4_1_pl = ht4_1_pl - u4_1 * u4_1_pl
                  h4_1_ph = ht4_1_ph - u4_1 * u4_1_ph
                  h4_1_mf = ht4_1_mf - u4_1 * u4_1_mf
                  h4_1_ml = ht4_1_ml - u4_1 * u4_1_ml
                  h4_1_mh = ht4_1_mh - u4_1 * u4_1_mh
                  h4_1_Tb = ht4_1_Tb - u4_1 * u4_1_Tb
                  h4_1_Mi = ht4_1_Mi - u4_1 * u4_1_Mi

                  Tguess = t4a + (h4_1 - h4a) / cp4a
                  T4_1, T4_1_h4_1, T_al = gas_tsetd(lambdap, nair, h4_1, Tguess)
                  T4_1_pl = T4_1_h4_1 * h4_1_pl
                  T4_1_ph = T4_1_h4_1 * h4_1_ph
                  T4_1_mf = T4_1_h4_1 * h4_1_mf
                  T4_1_ml = T4_1_h4_1 * h4_1_ml
                  T4_1_mh = T4_1_h4_1 * h4_1_mh
                  T4_1_Tb = T4_1_h4_1 * h4_1_Tb
                  T4_1_Mi = T4_1_h4_1 * h4_1_Mi
                  for i = 1:nair
                        T4_1_pl = T4_1_pl + T_al[i] * lamp_pl[i]
                        T4_1_ph = T4_1_ph + T_al[i] * lamp_ph[i]
                        T4_1_mf = T4_1_mf + T_al[i] * lamp_mf[i]
                        T4_1_ml = T4_1_ml + T_al[i] * lamp_ml[i]
                        T4_1_mh = T4_1_mh + T_al[i] * lamp_mh[i]
                        T4_1_Tb = T4_1_Tb + T_al[i] * lamp_Tb[i]
                        T4_1_Mi = T4_1_Mi + T_al[i] * lamp_Mi[i]
                  end
                  s4_1, s4_1_T4_1, h4_1, dum, cp4_1, R4_1 = gassum(lambdap, nair, T4_1)
                  s4_1_pl = s4_1_T4_1 * T4_1_pl
                  s4_1_ph = s4_1_T4_1 * T4_1_ph
                  s4_1_mf = s4_1_T4_1 * T4_1_mf
                  s4_1_ml = s4_1_T4_1 * T4_1_ml
                  s4_1_mh = s4_1_T4_1 * T4_1_mh
                  s4_1_Tb = s4_1_T4_1 * T4_1_Tb
                  s4_1_Mi = s4_1_T4_1 * T4_1_Mi
                  for i = 1:nair
                        si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, T4_1)
                        s4_1_pl = s4_1_pl + si * lamp_pl[i]
                        s4_1_ph = s4_1_ph + si * lamp_ph[i]
                        s4_1_mf = s4_1_mf + si * lamp_mf[i]
                        s4_1_ml = s4_1_ml + si * lamp_ml[i]
                        s4_1_mh = s4_1_mh + si * lamp_mh[i]
                        s4_1_Tb = s4_1_Tb + si * lamp_Tb[i]
                        s4_1_Mi = s4_1_Mi + si * lamp_Mi[i]
                  end

                  #----- all stagnation quantities, from total-static enthalpy difference
                  dhb = ht4_1 - h4_1
                  dhb_pl = ht4_1_pl - h4_1_pl
                  dhb_ph = ht4_1_ph - h4_1_ph
                  dhb_mf = ht4_1_mf - h4_1_mf
                  dhb_ml = ht4_1_ml - h4_1_ml
                  dhb_mh = ht4_1_mh - h4_1_mh
                  dhb_Tb = ht4_1_Tb - h4_1_Tb
                  dhb_Mi = ht4_1_Mi - h4_1_Mi
                  epi = 1.0
                  pt4_1, Tt4_1, ht4_1, st4_1, cpt4_1, Rt4_1,
                  pt4_1_s4_1,
                  pt4_1_p4_1,
                  dum,
                  pt4_1_h4_1, Tt4_1_h4_1, ht4_1_h4_1, st4_1_h4_1,
                  pt4_1_dhb, Tt4_1_dhb, ht4_1_dhb, st4_1_dhb,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_delhd(lambdap, nair,
                        p4_1, T4_1, h4_1, s4_1, cp4_1, R4_1, dhb, epi)

                  pt4_1_pl = pt4_1_s4_1 * s4_1_pl +
                            pt4_1_p4_1 * p4_1_pl +
                            pt4_1_h4_1 * h4_1_pl + pt4_1_dhb * dhb_pl

                  pt4_1_ph = pt4_1_s4_1 * s4_1_ph +
                            pt4_1_p4_1 * p4_1_ph +
                            pt4_1_h4_1 * h4_1_ph + pt4_1_dhb * dhb_ph

                  pt4_1_mf = pt4_1_s4_1 * s4_1_mf +
                            pt4_1_h4_1 * h4_1_mf + pt4_1_dhb * dhb_mf

                  pt4_1_ml = pt4_1_s4_1 * s4_1_ml +
                            pt4_1_p4_1 * p4_1_ml +
                            pt4_1_h4_1 * h4_1_ml + pt4_1_dhb * dhb_ml

                  pt4_1_mh = pt4_1_s4_1 * s4_1_mh +
                            pt4_1_p4_1 * p4_1_mh +
                            pt4_1_h4_1 * h4_1_mh + pt4_1_dhb * dhb_mh

                  pt4_1_Tb = pt4_1_s4_1 * s4_1_Tb +
                            pt4_1_p4_1 * p4_1_Tb +
                            pt4_1_h4_1 * h4_1_Tb + pt4_1_dhb * dhb_Tb

                  pt4_1_Mi = pt4_1_s4_1 * s4_1_Mi +
                            pt4_1_p4_1 * p4_1_Mi +
                            pt4_1_h4_1 * h4_1_Mi + pt4_1_dhb * dhb_Mi

                  for i = 1:nair
                        pt4_1_pl = pt4_1_pl + p_al[i] * lamp_pl[i]
                        pt4_1_ph = pt4_1_ph + p_al[i] * lamp_ph[i]
                        pt4_1_mf = pt4_1_mf + p_al[i] * lamp_mf[i]
                        pt4_1_ml = pt4_1_ml + p_al[i] * lamp_ml[i]
                        pt4_1_mh = pt4_1_mh + p_al[i] * lamp_mh[i]
                        pt4_1_Tb = pt4_1_Tb + p_al[i] * lamp_Tb[i]
                        pt4_1_Mi = pt4_1_Mi + p_al[i] * lamp_Mi[i]

                  end

            end

            # ===============================================================
            #---- HPT and LPT work

            #---- bypass ratio  (splitter = Splitter() constructed above Newton loop)
            BPR, BPR_mf_dir, BPR_ml_dir, BPR_pt2, BPR_pt1_9c =
                bypass_ratio(splitter, mf, ml, pt2, pt2ac, Tt2, Tt2ac)
            BPR_mf = BPR_mf_dir + BPR_pt2 * pt2_mf + BPR_pt1_9c * pt1_9c_mf
            BPR_ml = BPR_ml_dir + BPR_pt2 * pt2_ml + BPR_pt1_9c * pt1_9c_ml
            BPR_Mi =              BPR_pt2 * pt2_Mi  + BPR_pt1_9c * pt1_9c_Mi

            #---- HPT work with Jacobian  (shaft_hp = Shaft(epsh, 1.0) constructed above Newton loop)
            dhht, dhht_fo, dhht_ff, dhht_ht3, dhht_ht2_5c =
                hp_shaft_workd(shaft_hp, fo, ff, ht3, ht2_5c)

            dhht_pl = dhht_ff * ff_pl + dhht_ht3 * ht3_pl + dhht_ht2_5c * ht2_5c_pl
            dhht_ph = dhht_ff * ff_ph + dhht_ht3 * ht3_ph
            dhht_mf = dhht_ff * ff_mf + dhht_fo * fo_mf
            dhht_ml = dhht_ff * ff_ml + dhht_fo * fo_ml + dhht_ht3 * ht3_ml + dhht_ht2_5c * ht2_5c_ml
            dhht_mh = dhht_ff * ff_mh + dhht_ht3 * ht3_mh
            dhht_Tb = dhht_ff * ff_Tb
            dhht_Mi = dhht_ff * ff_Mi + dhht_fo * fo_Mi

            #---- LPT work with Jacobian  (shaft_lp = Shaft(epsl, Gearf) constructed above Newton loop)
            dhlt, dhlt_fo, dhlt_ff, dhlt_BPR,
            dhlt_ht2_5, dhlt_ht1_9c, dhlt_ht2_1, dhlt_ht2, dhlt_Pom,
            dlfac_fo_lp, dlfac_ff_lp =
                lp_shaft_workd(shaft_lp, fo, ff, BPR, ht2_5, ht2ac, ht13, ht2, Pom)

            dhlt_pf = dhlt_ht2_1 * ht2_1_pf
            dhlt_pl = dhlt_ff * ff_pl + dhlt_ht2_5 * ht2_5_pl
            dhlt_ph = dhlt_ff * ff_ph
            dhlt_mf = dhlt_ff * ff_mf + dhlt_fo * fo_mf +
                      dhlt_BPR * BPR_mf + dhlt_ht2_1 * ht2_1_mf + dhlt_Pom * Pom_mf
            # Mirrors upstream Fortran exactly: dhlt_ml uses demand_no_pom = demand − Pom
            # in the dlfac_ml contribution — a copy-paste omission in tfoper.f relative to
            # all other dhlt_* partials.  The correct partial (with +Pom) is tasopt-go7.
            dhlt_ml = (ht2_5 - ht2ac + BPR * (ht13 - ht2)) *
                          (dlfac_fo_lp * fo_ml + dlfac_ff_lp * ff_ml) +
                      dhlt_Pom * (ht2_5_ml + BPR_ml * (ht13 - ht2) + Pom_ml)
            dhlt_mh = dhlt_ff * ff_mh
            dhlt_Tb = dhlt_ff * ff_Tb
            dhlt_Mi = dhlt_ff * ff_Mi + dhlt_fo * fo_Mi +
                      dhlt_BPR * BPR_Mi + dhlt_Pom * Pom_Mi

            #---- HPT corrected mass flow, using LPC corrected mass flow and fuel/air ratio
            mbht = ml * (1.0 - fo + ff) * sqrt(Tt4_1 / Tt2ac) * pt2ac / pt4_1
            mbht_ml = (1.0 - fo + ff) * sqrt(Tt4_1 / Tt2ac) * pt2ac / pt4_1
            mbht_ff = ml * sqrt(Tt4_1 / Tt2ac) * pt2ac / pt4_1
            mbht_fo = -ml * sqrt(Tt4_1 / Tt2ac) * pt2ac / pt4_1
            mbht_Tt4_1 = 0.5 * mbht / Tt4_1
            mbht_pt4_1 = -mbht / pt4_1
            mbht_pt1_9c = ml * (1.0 - fo + ff) * sqrt(Tt4_1 / Tt2ac) / pt4_1

            mbht_pl = mbht_ff * ff_pl + mbht_Tt4_1 * Tt4_1_pl + mbht_pt4_1 * pt4_1_pl
            mbht_ph = mbht_ff * ff_ph + mbht_Tt4_1 * Tt4_1_ph + mbht_pt4_1 * pt4_1_ph
            mbht_mf = mbht_ff * ff_mf + mbht_Tt4_1 * Tt4_1_mf + mbht_pt4_1 * pt4_1_mf +
                      mbht_fo * fo_mf + mbht_pt1_9c * pt1_9c_mf
            mbht_ml = mbht_ff * ff_ml + mbht_Tt4_1 * Tt4_1_ml + mbht_pt4_1 * pt4_1_ml +
                      mbht_fo * fo_ml + mbht_pt1_9c * pt1_9c_ml + mbht_ml
            mbht_mh = mbht_ff * ff_mh + mbht_Tt4_1 * Tt4_1_mh + mbht_pt4_1 * pt4_1_mh
            mbht_Tb = mbht_ff * ff_Tb + mbht_Tt4_1 * Tt4_1_Tb + mbht_pt4_1 * pt4_1_Tb
            mbht_Mi = mbht_ff * ff_Mi + mbht_Tt4_1 * Tt4_1_Mi + mbht_pt4_1 * pt4_1_Mi +
                      mbht_fo * fo_Mi + mbht_pt1_9c * pt1_9c_Mi

            #---- HPT corrected speed
            Nbht = Nh * sqrt(Tt2_5c / Tt4_1)
            Nbht_Nh = sqrt(Tt2_5c / Tt4_1)
            Nbht_Tt2_5c = 0.5 * Nbht / Tt2_5c
            Nbht_Tt4_1 = -0.5 * Nbht / Tt4_1

            Nbht_pl = Nbht_Tt2_5c * Tt2_5c_pl + Nbht_Tt4_1 * Tt4_1_pl
            Nbht_ph = Nbht_Nh * Nh_ph + Nbht_Tt4_1 * Tt4_1_ph
            Nbht_mf = Nbht_Tt4_1 * Tt4_1_mf
            Nbht_ml = Nbht_Tt2_5c * Tt2_5c_ml + Nbht_Tt4_1 * Tt4_1_ml
            Nbht_mh = Nbht_Nh * Nh_mh + Nbht_Tt4_1 * Tt4_1_mh
            Nbht_Tb = Nbht_Tt4_1 * Tt4_1_Tb
            Nbht_Mi = Nbht_Tt4_1 * Tt4_1_Mi

            #---- HPT efficiency
            # HACK: HSC
            #First find uncooled HPT efficiency and derivatives
            epht1, epht1_dhht, epht1_mbht, epht1_Nbht, epht1_Tt4_1,
            epht1_cpt4_1, epht1_Rt4_1 = turbine_efficiency(turb_hp, dhht, mbht, Nbht, Tt4_1, cpt4_1, Rt4_1)

            #Find cooled HPT efficiency epht
            epht = find_cooled_hpt_efficiency(epht1, epht_fc, fc0, fc)

            if (epht < 0.80)
                  epht = 0.80
                  epht1_dhht = 0.0
                  epht1_mbht = 0.0
                  epht1_Nbht = 0.0
                  epht1_Tt4_1 = 0.0
                  epht1_cpt4_1 = 0.0
                  epht1_Rt4_1 = 0.0
            end

            epht1_pl = epht1_dhht * dhht_pl + epht1_mbht * mbht_pl + epht1_Nbht * Nbht_pl
            epht1_ph = epht1_dhht * dhht_ph + epht1_mbht * mbht_ph + epht1_Nbht * Nbht_ph
            epht1_mf = epht1_dhht * dhht_mf + epht1_mbht * mbht_mf + epht1_Nbht * Nbht_mf
            epht1_ml = epht1_dhht * dhht_ml + epht1_mbht * mbht_ml + epht1_Nbht * Nbht_ml
            epht1_mh = epht1_dhht * dhht_mh + epht1_mbht * mbht_mh + epht1_Nbht * Nbht_mh
            epht1_Tb = epht1_dhht * dhht_Tb + epht1_mbht * mbht_Tb + epht1_Nbht * Nbht_Tb
            epht1_Mi = epht1_dhht * dhht_Mi + epht1_mbht * mbht_Mi + epht1_Nbht * Nbht_Mi

            epht1_pl = epht1_Tt4_1 * Tt4_1_pl + epht1_cpt4_1 * cpt4_1_pl +
                      epht1_Rt4_1 * Rt4_1_pl + epht1_pl
            epht1_ph = epht1_Tt4_1 * Tt4_1_ph + epht1_cpt4_1 * cpt4_1_ph +
                      epht1_Rt4_1 * Rt4_1_ph + epht1_ph
            epht1_mf = epht1_Tt4_1 * Tt4_1_mf + epht1_cpt4_1 * cpt4_1_mf +
                      epht1_Rt4_1 * Rt4_1_mf + epht1_mf
            epht1_ml = epht1_Tt4_1 * Tt4_1_ml + epht1_cpt4_1 * cpt4_1_ml +
                      epht1_Rt4_1 * Rt4_1_ml + epht1_ml
            epht1_mh = epht1_Tt4_1 * Tt4_1_mh + epht1_cpt4_1 * cpt4_1_mh +
                      epht1_Rt4_1 * Rt4_1_mh + epht1_mh
            epht1_Tb = epht1_Tt4_1 * Tt4_1_Tb + epht1_cpt4_1 * cpt4_1_Tb +
                      epht1_Rt4_1 * Rt4_1_Tb + epht1_Tb
            epht1_Mi = epht1_Tt4_1 * Tt4_1_Mi + epht1_cpt4_1 * cpt4_1_Mi +
                      epht1_Rt4_1 * Rt4_1_Mi + epht1_Mi

            #Chain rule for HPT efficiency derivatives
            epht_pl = epht1_pl + epht_fc * fc_pl
            epht_ph = epht1_ph + epht_fc * fc_ph
            epht_mf = epht1_mf + epht_fc * fc_mf
            epht_ml = epht1_ml + epht_fc * fc_ml
            epht_mh = epht1_mh + epht_fc * fc_mh
            epht_Tb = epht1_Tb + epht_fc * fc_Tb
            epht_Mi = epht1_Mi + epht_fc * fc_Mi
           
            #---- HPT work to determine station 4.5
            pt4_5, Tt4_5, ht4_5, st4_5, cpt4_5, Rt4_5,
            pt4_5_pt4_1, pt4_5_st4_1,
            pt4_5_ht4_1, Tt4_5_ht4_1, ht4_5_ht4_1, st4_5_ht4_1,
            pt4_5_dhht, Tt4_5_dhht, ht4_5_dhht, st4_5_dhht,
            pt4_5_epht,
            p_al, T_al, h_al, s_al = turbine_delhd(lambdap, nair,
                  pt4_1, Tt4_1, ht4_1, st4_1, cpt4_1, Rt4_1, dhht, epht)

            pt4_5_pl = pt4_5_st4_1 * st4_1_pl +
                      pt4_5_pt4_1 * pt4_1_pl +
                      pt4_5_epht * epht_pl +
                      pt4_5_ht4_1 * ht4_1_pl + pt4_5_dhht * dhht_pl
            Tt4_5_pl = Tt4_5_ht4_1 * ht4_1_pl + Tt4_5_dhht * dhht_pl
            ht4_5_pl = ht4_5_ht4_1 * ht4_1_pl + ht4_5_dhht * dhht_pl
            st4_5_pl = st4_5_ht4_1 * ht4_1_pl + st4_5_dhht * dhht_pl

            pt4_5_ph = pt4_5_st4_1 * st4_1_ph +
                      pt4_5_pt4_1 * pt4_1_ph +
                      pt4_5_epht * epht_ph +
                      pt4_5_ht4_1 * ht4_1_ph + pt4_5_dhht * dhht_ph
            Tt4_5_ph = Tt4_5_ht4_1 * ht4_1_ph + Tt4_5_dhht * dhht_ph
            ht4_5_ph = ht4_5_ht4_1 * ht4_1_ph + ht4_5_dhht * dhht_ph
            st4_5_ph = st4_5_ht4_1 * ht4_1_ph + st4_5_dhht * dhht_ph

            pt4_5_mf = pt4_5_st4_1 * st4_1_mf +
                      pt4_5_pt4_1 * pt4_1_mf +
                      pt4_5_epht * epht_mf +
                      pt4_5_ht4_1 * ht4_1_mf + pt4_5_dhht * dhht_mf
            Tt4_5_mf = Tt4_5_ht4_1 * ht4_1_mf + Tt4_5_dhht * dhht_mf
            ht4_5_mf = ht4_5_ht4_1 * ht4_1_mf + ht4_5_dhht * dhht_mf
            st4_5_mf = st4_5_ht4_1 * ht4_1_mf + st4_5_dhht * dhht_mf

            pt4_5_ml = pt4_5_st4_1 * st4_1_ml +
                      pt4_5_pt4_1 * pt4_1_ml +
                      pt4_5_epht * epht_ml +
                      pt4_5_ht4_1 * ht4_1_ml + pt4_5_dhht * dhht_ml
            Tt4_5_ml = Tt4_5_ht4_1 * ht4_1_ml + Tt4_5_dhht * dhht_ml
            ht4_5_ml = ht4_5_ht4_1 * ht4_1_ml + ht4_5_dhht * dhht_ml
            st4_5_ml = st4_5_ht4_1 * ht4_1_ml + st4_5_dhht * dhht_ml

            pt4_5_mh = pt4_5_st4_1 * st4_1_mh +
                      pt4_5_pt4_1 * pt4_1_mh +
                      pt4_5_epht * epht_mh +
                      pt4_5_ht4_1 * ht4_1_mh + pt4_5_dhht * dhht_mh
            Tt4_5_mh = Tt4_5_ht4_1 * ht4_1_mh + Tt4_5_dhht * dhht_mh
            ht4_5_mh = ht4_5_ht4_1 * ht4_1_mh + ht4_5_dhht * dhht_mh
            st4_5_mh = st4_5_ht4_1 * ht4_1_mh + st4_5_dhht * dhht_mh

            pt4_5_Tb = pt4_5_st4_1 * st4_1_Tb +
                      pt4_5_pt4_1 * pt4_1_Tb +
                      pt4_5_epht * epht_Tb +
                      pt4_5_ht4_1 * ht4_1_Tb + pt4_5_dhht * dhht_Tb
            Tt4_5_Tb = Tt4_5_ht4_1 * ht4_1_Tb + Tt4_5_dhht * dhht_Tb
            ht4_5_Tb = ht4_5_ht4_1 * ht4_1_Tb + ht4_5_dhht * dhht_Tb
            st4_5_Tb = st4_5_ht4_1 * ht4_1_Tb + st4_5_dhht * dhht_Tb

            pt4_5_Mi = pt4_5_st4_1 * st4_1_Mi +
                      pt4_5_pt4_1 * pt4_1_Mi +
                      pt4_5_epht * epht_Mi +
                      pt4_5_ht4_1 * ht4_1_Mi + pt4_5_dhht * dhht_Mi
            Tt4_5_Mi = Tt4_5_ht4_1 * ht4_1_Mi + Tt4_5_dhht * dhht_Mi
            ht4_5_Mi = ht4_5_ht4_1 * ht4_1_Mi + ht4_5_dhht * dhht_Mi
            st4_5_Mi = st4_5_ht4_1 * ht4_1_Mi + st4_5_dhht * dhht_Mi


            for i = 1:nair
                  pt4_5_pl = pt4_5_pl + p_al[i] * lamp_pl[i]
                  Tt4_5_pl = Tt4_5_pl + T_al[i] * lamp_pl[i]
                  ht4_5_pl = ht4_5_pl + h_al[i] * lamp_pl[i]
                  st4_5_pl = st4_5_pl + s_al[i] * lamp_pl[i]

                  pt4_5_ph = pt4_5_ph + p_al[i] * lamp_ph[i]
                  Tt4_5_ph = Tt4_5_ph + T_al[i] * lamp_ph[i]
                  ht4_5_ph = ht4_5_ph + h_al[i] * lamp_ph[i]
                  st4_5_ph = st4_5_ph + s_al[i] * lamp_ph[i]

                  pt4_5_mf = pt4_5_mf + p_al[i] * lamp_mf[i]
                  Tt4_5_mf = Tt4_5_mf + T_al[i] * lamp_mf[i]
                  ht4_5_mf = ht4_5_mf + h_al[i] * lamp_mf[i]
                  st4_5_mf = st4_5_mf + s_al[i] * lamp_mf[i]

                  pt4_5_ml = pt4_5_ml + p_al[i] * lamp_ml[i]
                  Tt4_5_ml = Tt4_5_ml + T_al[i] * lamp_ml[i]
                  ht4_5_ml = ht4_5_ml + h_al[i] * lamp_ml[i]
                  st4_5_ml = st4_5_ml + s_al[i] * lamp_ml[i]

                  pt4_5_mh = pt4_5_mh + p_al[i] * lamp_mh[i]
                  Tt4_5_mh = Tt4_5_mh + T_al[i] * lamp_mh[i]
                  ht4_5_mh = ht4_5_mh + h_al[i] * lamp_mh[i]
                  st4_5_mh = st4_5_mh + s_al[i] * lamp_mh[i]

                  pt4_5_Tb = pt4_5_Tb + p_al[i] * lamp_Tb[i]
                  Tt4_5_Tb = Tt4_5_Tb + T_al[i] * lamp_Tb[i]
                  ht4_5_Tb = ht4_5_Tb + h_al[i] * lamp_Tb[i]
                  st4_5_Tb = st4_5_Tb + s_al[i] * lamp_Tb[i]

                  pt4_5_Mi = pt4_5_Mi + p_al[i] * lamp_Mi[i]
                  Tt4_5_Mi = Tt4_5_Mi + T_al[i] * lamp_Mi[i]
                  ht4_5_Mi = ht4_5_Mi + h_al[i] * lamp_Mi[i]
                  st4_5_Mi = st4_5_Mi + s_al[i] * lamp_Mi[i]
            end

            #---- will also need cpt4_5,Rt4_5 derivatives
            st4_5, st4_5_Tt4_5,
            ht4_5, ht4_5_Tt4_5,
            cpt4_5, cpt4_5_Tt4_5, Rt4_5 = gassumd(lambdap, nair, Tt4_5)
            cpt4_5_pl = cpt4_5_Tt4_5 * Tt4_5_pl
            cpt4_5_ph = cpt4_5_Tt4_5 * Tt4_5_ph
            cpt4_5_mf = cpt4_5_Tt4_5 * Tt4_5_mf
            cpt4_5_ml = cpt4_5_Tt4_5 * Tt4_5_ml
            cpt4_5_mh = cpt4_5_Tt4_5 * Tt4_5_mh
            cpt4_5_Tb = cpt4_5_Tt4_5 * Tt4_5_Tb
            cpt4_5_Mi = cpt4_5_Tt4_5 * Tt4_5_Mi

            Rt4_5_pl = 0.0
            Rt4_5_ph = 0.0
            Rt4_5_mf = 0.0
            Rt4_5_ml = 0.0
            Rt4_5_mh = 0.0
            Rt4_5_Tb = 0.0
            Rt4_5_Mi = 0.0
            for i = 1:nair
                  si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, Tt4_5)
                  cpt4_5_pl = cpt4_5_pl + cpi * lamp_pl[i]
                  cpt4_5_ph = cpt4_5_ph + cpi * lamp_ph[i]
                  cpt4_5_mf = cpt4_5_mf + cpi * lamp_mf[i]
                  cpt4_5_ml = cpt4_5_ml + cpi * lamp_ml[i]
                  cpt4_5_mh = cpt4_5_mh + cpi * lamp_mh[i]
                  cpt4_5_Tb = cpt4_5_Tb + cpi * lamp_Tb[i]
                  cpt4_5_Mi = cpt4_5_Mi + cpi * lamp_Mi[i]

                  Rt4_5_pl = Rt4_5_pl + Ri * lamp_pl[i]
                  Rt4_5_ph = Rt4_5_ph + Ri * lamp_ph[i]
                  Rt4_5_mf = Rt4_5_mf + Ri * lamp_mf[i]
                  Rt4_5_ml = Rt4_5_ml + Ri * lamp_ml[i]
                  Rt4_5_mh = Rt4_5_mh + Ri * lamp_mh[i]
                  Rt4_5_Tb = Rt4_5_Tb + Ri * lamp_Tb[i]
                  Rt4_5_Mi = Rt4_5_Mi + Ri * lamp_Mi[i]
            end


            #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            #---- LPT corrected mass flow (in terms of Newton variables)
            mblt = ml * (1.0 - fo + ff) * sqrt(Tt4_5 / Tt2ac) * pt2ac / pt4_5
            mblt_ml = (1.0 - fo + ff) * sqrt(Tt4_5 / Tt2ac) * pt2ac / pt4_5
            mblt_ff = ml * sqrt(Tt4_5 / Tt2ac) * pt2ac / pt4_5
            mblt_fo = -ml * sqrt(Tt4_5 / Tt2ac) * pt2ac / pt4_5
            mblt_Tt4_5 = 0.5 * mblt / Tt4_5
            mblt_pt4_5 = -mblt / pt4_5
            mblt_pt1_9c = ml * (1.0 - fo + ff) * sqrt(Tt4_5 / Tt2ac) / pt4_5

            mblt_pl = mblt_ff * ff_pl + mblt_Tt4_5 * Tt4_5_pl + mblt_pt4_5 * pt4_5_pl
            mblt_ph = mblt_ff * ff_ph + mblt_Tt4_5 * Tt4_5_ph + mblt_pt4_5 * pt4_5_ph
            mblt_mf = mblt_ff * ff_mf + mblt_Tt4_5 * Tt4_5_mf + mblt_pt4_5 * pt4_5_mf +
                      mblt_fo * fo_mf + mblt_pt1_9c * pt1_9c_mf
            mblt_ml = mblt_ff * ff_ml + mblt_Tt4_5 * Tt4_5_ml + mblt_pt4_5 * pt4_5_ml +
                      mblt_fo * fo_ml + mblt_pt1_9c * pt1_9c_ml + mblt_ml
            mblt_mh = mblt_ff * ff_mh + mblt_Tt4_5 * Tt4_5_mh + mblt_pt4_5 * pt4_5_mh
            mblt_Tb = mblt_ff * ff_Tb + mblt_Tt4_5 * Tt4_5_Tb + mblt_pt4_5 * pt4_5_Tb
            mblt_Mi = mblt_ff * ff_Mi + mblt_Tt4_5 * Tt4_5_Mi + mblt_pt4_5 * pt4_5_Mi +
                      mblt_fo * fo_Mi + mblt_pt1_9c * pt1_9c_Mi

            #---- LPT corrected speed
            Nblt = Nl * sqrt(Tt2ac / Tt4_5)
            Nblt_Nl = sqrt(Tt2ac / Tt4_5)
            Nblt_Tt4_5 = -0.5 * Nblt / Tt4_5

            Nblt_pl = Nblt_Nl * Nl_pl + Nblt_Tt4_5 * Tt4_5_pl
            Nblt_ph = Nblt_Tt4_5 * Tt4_5_ph
            Nblt_mf = Nblt_Tt4_5 * Tt4_5_mf
            Nblt_ml = Nblt_Nl * Nl_ml + Nblt_Tt4_5 * Tt4_5_ml
            Nblt_mh = Nblt_Tt4_5 * Tt4_5_mh
            Nblt_Tb = Nblt_Tt4_5 * Tt4_5_Tb
            Nblt_Mi = Nblt_Tt4_5 * Tt4_5_Mi

            #---- LPT efficiency
            eplt,
            eplt_dhlt, eplt_mblt, eplt_Nblt,
            eplt_Tt4_5, eplt_cpt4_5, eplt_Rt4_5 = turbine_efficiency(turb_lp, dhlt, mblt, Nblt, Tt4_5, cpt4_5, Rt4_5)

            if (eplt < 0.80)
                  eplt = 0.80
                  eplt_dhlt = 0.0
                  eplt_mblt = 0.0
                  eplt_Nblt = 0.0
                  eplt_Tt4_5 = 0.0
                  eplt_cpt4_5 = 0.0
                  eplt_Rt4_5 = 0.0
            end

            eplt_pf = eplt_dhlt * dhlt_pf
            eplt_pl = eplt_dhlt * dhlt_pl + eplt_mblt * mblt_pl + eplt_Nblt * Nblt_pl
            eplt_ph = eplt_dhlt * dhlt_ph + eplt_mblt * mblt_ph + eplt_Nblt * Nblt_ph
            eplt_mf = eplt_dhlt * dhlt_mf + eplt_mblt * mblt_mf + eplt_Nblt * Nblt_mf
            eplt_ml = eplt_dhlt * dhlt_ml + eplt_mblt * mblt_ml + eplt_Nblt * Nblt_ml
            eplt_mh = eplt_dhlt * dhlt_mh + eplt_mblt * mblt_mh + eplt_Nblt * Nblt_mh
            eplt_Tb = eplt_dhlt * dhlt_Tb + eplt_mblt * mblt_Tb + eplt_Nblt * Nblt_Tb
            eplt_Mi = eplt_dhlt * dhlt_Mi + eplt_mblt * mblt_Mi + eplt_Nblt * Nblt_Mi

            eplt_pl = eplt_Tt4_5 * Tt4_5_pl + eplt_cpt4_5 * cpt4_5_pl +
                      eplt_Rt4_5 * Rt4_5_pl + eplt_pl
            eplt_ph = eplt_Tt4_5 * Tt4_5_ph + eplt_cpt4_5 * cpt4_5_ph +
                      eplt_Rt4_5 * Rt4_5_ph + eplt_ph
            eplt_mf = eplt_Tt4_5 * Tt4_5_mf + eplt_cpt4_5 * cpt4_5_mf +
                      eplt_Rt4_5 * Rt4_5_mf + eplt_mf
            eplt_ml = eplt_Tt4_5 * Tt4_5_ml + eplt_cpt4_5 * cpt4_5_ml +
                      eplt_Rt4_5 * Rt4_5_ml + eplt_ml
            eplt_mh = eplt_Tt4_5 * Tt4_5_mh + eplt_cpt4_5 * cpt4_5_mh +
                      eplt_Rt4_5 * Rt4_5_mh + eplt_mh
            eplt_Tb = eplt_Tt4_5 * Tt4_5_Tb + eplt_cpt4_5 * cpt4_5_Tb +
                      eplt_Rt4_5 * Rt4_5_Tb + eplt_Tb
            eplt_Mi = eplt_Tt4_5 * Tt4_5_Mi + eplt_cpt4_5 * cpt4_5_Mi +
                      eplt_Rt4_5 * Rt4_5_Mi + eplt_Mi


            if (Pc == 0.0)
                  #----- initial guess for pt8
                  epi = 1.0 / eplt
                  pt5, Tt5, ht5, st5, cpt5, Rt5 = gas_delh(lambdap, nair,
                        pt4_5, Tt4_5, ht4_5, st4_5, cpt4_5, Rt4_5, dhlt, epi)
                  Pc = pt5 * pitn
                  Pc = max(Pc, p0 * (1.0 + 0.2 * M0^2)^3.5)
            end

            Pc = max(Pc, 1.000001 * p0)

            pilt = Pc / pt4_5 / pitn
            pilt_Pc = 1.0 / pt4_5 / pitn
            pilt_pt4_5 = -pilt / pt4_5

            epi = 1.0 / eplt
            epi_eplt = -epi / eplt
            pt5, Tt5, ht5, st5, cpt5, Rt5,
            pt4_9_pt4_5,
            pt4_9_st4_5, Tt4_9_st4_5, ht4_9_st4_5, st4_9_st4_5,
            pt4_9_pilt, Tt4_9_pilt, ht4_9_pilt, st4_9_pilt,
            pt4_9_epi, Tt4_9_epi, ht4_9_epi, st4_9_epi,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_pratd(lambdap, nair,
                  pt4_5, Tt4_5, ht4_5, st4_5, cpt4_5, Rt4_5, pilt, epi)

            pt4_9_Pc = pt4_9_pilt * pilt_Pc
            Tt4_9_Pc = Tt4_9_pilt * pilt_Pc
            ht4_9_Pc = ht4_9_pilt * pilt_Pc
            st4_9_Pc = st4_9_pilt * pilt_Pc

            pt4_9_pt4_5 = pt4_9_pilt * pilt_pt4_5 + pt4_9_pt4_5
            Tt4_9_pt4_5 = Tt4_9_pilt * pilt_pt4_5
            ht4_9_pt4_5 = ht4_9_pilt * pilt_pt4_5
            st4_9_pt4_5 = st4_9_pilt * pilt_pt4_5

            pt4_9_eplt = pt4_9_epi * epi_eplt
            Tt4_9_eplt = Tt4_9_epi * epi_eplt
            ht4_9_eplt = ht4_9_epi * epi_eplt
            st4_9_eplt = st4_9_epi * epi_eplt


            pt4_9_pf = pt4_9_eplt * eplt_pf
            Tt4_9_pf = Tt4_9_eplt * eplt_pf
            ht4_9_pf = ht4_9_eplt * eplt_pf
            st4_9_pf = st4_9_eplt * eplt_pf

            pt4_9_pl = pt4_9_pt4_5 * pt4_5_pl + pt4_9_st4_5 * st4_5_pl + pt4_9_eplt * eplt_pl
            Tt4_9_pl = Tt4_9_pt4_5 * pt4_5_pl + Tt4_9_st4_5 * st4_5_pl + Tt4_9_eplt * eplt_pl
            ht4_9_pl = ht4_9_pt4_5 * pt4_5_pl + ht4_9_st4_5 * st4_5_pl + ht4_9_eplt * eplt_pl
            st4_9_pl = st4_9_pt4_5 * pt4_5_pl + st4_9_st4_5 * st4_5_pl + st4_9_eplt * eplt_pl

            pt4_9_ph = pt4_9_pt4_5 * pt4_5_ph + pt4_9_st4_5 * st4_5_ph + pt4_9_eplt * eplt_ph
            Tt4_9_ph = Tt4_9_pt4_5 * pt4_5_ph + Tt4_9_st4_5 * st4_5_ph + Tt4_9_eplt * eplt_ph
            ht4_9_ph = ht4_9_pt4_5 * pt4_5_ph + ht4_9_st4_5 * st4_5_ph + ht4_9_eplt * eplt_ph
            st4_9_ph = st4_9_pt4_5 * pt4_5_ph + st4_9_st4_5 * st4_5_ph + st4_9_eplt * eplt_ph

            pt4_9_mf = pt4_9_pt4_5 * pt4_5_mf + pt4_9_st4_5 * st4_5_mf + pt4_9_eplt * eplt_mf
            Tt4_9_mf = Tt4_9_pt4_5 * pt4_5_mf + Tt4_9_st4_5 * st4_5_mf + Tt4_9_eplt * eplt_mf
            ht4_9_mf = ht4_9_pt4_5 * pt4_5_mf + ht4_9_st4_5 * st4_5_mf + ht4_9_eplt * eplt_mf
            st4_9_mf = st4_9_pt4_5 * pt4_5_mf + st4_9_st4_5 * st4_5_mf + st4_9_eplt * eplt_mf

            pt4_9_ml = pt4_9_pt4_5 * pt4_5_ml + pt4_9_st4_5 * st4_5_ml + pt4_9_eplt * eplt_ml
            Tt4_9_ml = Tt4_9_pt4_5 * pt4_5_ml + Tt4_9_st4_5 * st4_5_ml + Tt4_9_eplt * eplt_ml
            ht4_9_ml = ht4_9_pt4_5 * pt4_5_ml + ht4_9_st4_5 * st4_5_ml + ht4_9_eplt * eplt_ml
            st4_9_ml = st4_9_pt4_5 * pt4_5_ml + st4_9_st4_5 * st4_5_ml + st4_9_eplt * eplt_ml

            pt4_9_mh = pt4_9_pt4_5 * pt4_5_mh + pt4_9_st4_5 * st4_5_mh + pt4_9_eplt * eplt_mh
            Tt4_9_mh = Tt4_9_pt4_5 * pt4_5_mh + Tt4_9_st4_5 * st4_5_mh + Tt4_9_eplt * eplt_mh
            ht4_9_mh = ht4_9_pt4_5 * pt4_5_mh + ht4_9_st4_5 * st4_5_mh + ht4_9_eplt * eplt_mh
            st4_9_mh = st4_9_pt4_5 * pt4_5_mh + st4_9_st4_5 * st4_5_mh + st4_9_eplt * eplt_mh


            pt4_9_Tb = pt4_9_pt4_5 * pt4_5_Tb + pt4_9_st4_5 * st4_5_Tb + pt4_9_eplt * eplt_Tb
            Tt4_9_Tb = Tt4_9_pt4_5 * pt4_5_Tb + Tt4_9_st4_5 * st4_5_Tb + Tt4_9_eplt * eplt_Tb
            ht4_9_Tb = ht4_9_pt4_5 * pt4_5_Tb + ht4_9_st4_5 * st4_5_Tb + ht4_9_eplt * eplt_Tb
            st4_9_Tb = st4_9_pt4_5 * pt4_5_Tb + st4_9_st4_5 * st4_5_Tb + st4_9_eplt * eplt_Tb

            pt4_9_Mi = pt4_9_pt4_5 * pt4_5_Mi + pt4_9_st4_5 * st4_5_Mi + pt4_9_eplt * eplt_Mi
            Tt4_9_Mi = Tt4_9_pt4_5 * pt4_5_Mi + Tt4_9_st4_5 * st4_5_Mi + Tt4_9_eplt * eplt_Mi
            ht4_9_Mi = ht4_9_pt4_5 * pt4_5_Mi + ht4_9_st4_5 * st4_5_Mi + ht4_9_eplt * eplt_Mi
            st4_9_Mi = st4_9_pt4_5 * pt4_5_Mi + st4_9_st4_5 * st4_5_Mi + st4_9_eplt * eplt_Mi

            for i = 1:nair
                  pt4_9_pl = pt4_9_pl + p_al[i] * lamp_pl[i]
                  Tt4_9_pl = Tt4_9_pl + T_al[i] * lamp_pl[i]
                  ht4_9_pl = ht4_9_pl + h_al[i] * lamp_pl[i]
                  st4_9_pl = st4_9_pl + s_al[i] * lamp_pl[i]

                  pt4_9_ph = pt4_9_ph + p_al[i] * lamp_ph[i]
                  Tt4_9_ph = Tt4_9_ph + T_al[i] * lamp_ph[i]
                  ht4_9_ph = ht4_9_ph + h_al[i] * lamp_ph[i]
                  st4_9_ph = st4_9_ph + s_al[i] * lamp_ph[i]

                  pt4_9_mf = pt4_9_mf + p_al[i] * lamp_mf[i]
                  Tt4_9_mf = Tt4_9_mf + T_al[i] * lamp_mf[i]
                  ht4_9_mf = ht4_9_mf + h_al[i] * lamp_mf[i]
                  st4_9_mf = st4_9_mf + s_al[i] * lamp_mf[i]

                  pt4_9_ml = pt4_9_ml + p_al[i] * lamp_ml[i]
                  Tt4_9_ml = Tt4_9_ml + T_al[i] * lamp_ml[i]
                  ht4_9_ml = ht4_9_ml + h_al[i] * lamp_ml[i]
                  st4_9_ml = st4_9_ml + s_al[i] * lamp_ml[i]

                  pt4_9_mh = pt4_9_mh + p_al[i] * lamp_mh[i]
                  Tt4_9_mh = Tt4_9_mh + T_al[i] * lamp_mh[i]
                  ht4_9_mh = ht4_9_mh + h_al[i] * lamp_mh[i]
                  st4_9_mh = st4_9_mh + s_al[i] * lamp_mh[i]

                  pt4_9_Tb = pt4_9_Tb + p_al[i] * lamp_Tb[i]
                  Tt4_9_Tb = Tt4_9_Tb + T_al[i] * lamp_Tb[i]
                  ht4_9_Tb = ht4_9_Tb + h_al[i] * lamp_Tb[i]
                  st4_9_Tb = st4_9_Tb + s_al[i] * lamp_Tb[i]

                  pt4_9_Mi = pt4_9_Mi + p_al[i] * lamp_Mi[i]
                  Tt4_9_Mi = Tt4_9_Mi + T_al[i] * lamp_Mi[i]
                  ht4_9_Mi = ht4_9_Mi + h_al[i] * lamp_Mi[i]
                  st4_9_Mi = st4_9_Mi + s_al[i] * lamp_Mi[i]
            end
            Tt4_9_ht4_9 = 1 / cpt5
            st4_9_Tt4_9 = cpt5 / Tt5

            # ===============================================================
            #---- Regenerative cooling heat exchanger 49-49c
            pt5c = pt5 - Δp_Regen
            ht5c = ht5 + Δh_Regen
            Tt5c, Tt4_9c_ht4_9c, _ = gas_tsetd(lambdap, nair, ht5c, Tt5)
            st5c, st4_9c_Tt4_9c, ht5c, ht4_9c_Tt4_9, cpt5c, cpt4_9c_Tt4_9c, Rt5c = gassumd(lambdap, nair, Tt5c)

            #Derivatives with respect to pressure are unchanged
            pt4_9c_pt4_9 = 1.0
            ht4_9c_ht4_9 = 1.0

            Tt4_9c_Tt4_9 = Tt4_9c_ht4_9c * ht4_9c_ht4_9 / Tt4_9_ht4_9
            st4_9c_st4_9 = st4_9c_Tt4_9c * Tt4_9c_Tt4_9 / st4_9_Tt4_9

            #Apply chain rule
            pt4_9c_pf = pt4_9_pf * pt4_9c_pt4_9
            Tt4_9c_pf = Tt4_9_pf * pt4_9c_pt4_9
            ht4_9c_pf = ht4_9_pf * pt4_9c_pt4_9
            st4_9c_pf = st4_9_pf * pt4_9c_pt4_9

            pt4_9c_pl = pt4_9_pl * pt4_9c_pt4_9 
            Tt4_9c_pl = Tt4_9_pl * Tt4_9c_Tt4_9
            ht4_9c_pl = ht4_9_pl * ht4_9c_ht4_9
            st4_9c_pl = st4_9_pl * st4_9c_st4_9

            pt4_9c_ph = pt4_9_ph * pt4_9c_pt4_9  
            Tt4_9c_ph = Tt4_9_ph * Tt4_9c_Tt4_9
            ht4_9c_ph = ht4_9_ph * ht4_9c_ht4_9
            st4_9c_ph = st4_9_ph * st4_9c_st4_9

            pt4_9c_mf = pt4_9_mf * pt4_9c_pt4_9  
            Tt4_9c_mf = Tt4_9_mf * Tt4_9c_Tt4_9
            ht4_9c_mf = ht4_9_mf * ht4_9c_ht4_9
            st4_9c_mf = st4_9_mf * st4_9c_st4_9

            pt4_9c_ml = pt4_9_ml * pt4_9c_pt4_9 
            Tt4_9c_ml = Tt4_9_ml * Tt4_9c_Tt4_9
            ht4_9c_ml = ht4_9_ml * ht4_9c_ht4_9
            st4_9c_ml = st4_9_ml * st4_9c_st4_9

            pt4_9c_mh = pt4_9_mh * pt4_9c_pt4_9  
            Tt4_9c_mh = Tt4_9_mh * Tt4_9c_Tt4_9
            ht4_9c_mh = ht4_9_mh * ht4_9c_ht4_9
            st4_9c_mh = st4_9_mh * st4_9c_st4_9

            pt4_9c_Tb = pt4_9_Tb * pt4_9c_pt4_9  
            Tt4_9c_Tb = Tt4_9_Tb * Tt4_9c_Tt4_9
            ht4_9c_Tb = ht4_9_Tb * ht4_9c_ht4_9
            st4_9c_Tb = st4_9_Tb * st4_9c_st4_9

            pt4_9c_Pc = pt4_9_Pc * pt4_9c_pt4_9  
            Tt4_9c_Pc = Tt4_9_Pc * Tt4_9c_Tt4_9
            ht4_9c_Pc = ht4_9_Pc * ht4_9c_ht4_9
            st4_9c_Pc = st4_9_Pc * st4_9c_st4_9

            pt4_9c_Mi = pt4_9_Mi * pt4_9c_pt4_9  
            Tt4_9c_Mi = Tt4_9_Mi * Tt4_9c_Tt4_9
            ht4_9c_Mi = ht4_9_Mi * ht4_9c_ht4_9
            st4_9c_Mi = st4_9_Mi * st4_9c_st4_9

            #---- apply core nozzle loss
            pt8 = pt5c * pitn
            Tt8 = Tt5c
            ht8 = ht5c
            st8 = st5c
            cpt8 = cpt5c
            Rt8 = Rt5c

            pt5_pf = pt4_9c_pf * pitn
            Tt5_pf = Tt4_9c_pf
            ht5_pf = ht4_9c_pf
            st5_pf = st4_9c_pf

            pt5_pl = pt4_9c_pl * pitn
            Tt5_pl = Tt4_9c_pl
            ht5_pl = ht4_9c_pl
            st5_pl = st4_9c_pl

            pt5_ph = pt4_9c_ph * pitn
            Tt5_ph = Tt4_9c_ph
            ht5_ph = ht4_9c_ph
            st5_ph = st4_9c_ph

            pt5_mf = pt4_9c_mf * pitn
            Tt5_mf = Tt4_9c_mf
            ht5_mf = ht4_9c_mf
            st5_mf = st4_9c_mf

            pt5_ml = pt4_9c_ml * pitn
            Tt5_ml = Tt4_9c_ml
            ht5_ml = ht4_9c_ml
            st5_ml = st4_9c_ml

            pt5_mh = pt4_9c_mh * pitn
            Tt5_mh = Tt4_9c_mh
            ht5_mh = ht4_9c_mh
            st5_mh = st4_9c_mh

            pt5_Tb = pt4_9c_Tb * pitn
            Tt5_Tb = Tt4_9c_Tb
            ht5_Tb = ht4_9c_Tb
            st5_Tb = st4_9c_Tb

            pt5_Pc = pt4_9c_Pc * pitn
            Tt5_Pc = Tt4_9c_Pc
            ht5_Pc = ht4_9c_Pc
            st5_Pc = st4_9c_Pc

            pt5_Mi = pt4_9c_Mi * pitn
            Tt5_Mi = Tt4_9c_Mi
            ht5_Mi = ht4_9c_Mi
            st5_Mi = st4_9c_Mi

            # ===============================================================
            #---- fan nozzle exit state (for function outputs and downstream thrust)
            #     nozzle_fan = Nozzle(one(T), A7) constructed above Newton loop
            p7, T7, h7, s7, cp7, R7, u7, rho7, M7 =
                  nozzle_exit(nozzle_fan, alpha, nair, pt18, Tt18, ht18, st18, cpt18, Rt18, p0)

            # ===============================================================
            #---- core nozzle exit state (for function outputs and species chain-rule)
            #     nozzle_core = Nozzle(one(T), A5) constructed above Newton loop
            p5, T5, h5, s5, cp5, R5, u5, rho5, M5 =
                  nozzle_exit(nozzle_core, lambdap, nair, pt8, Tt8, ht8, st8, cpt8, Rt8, p0)

            # ===========================================================================
            #---- set up Newton system

            #---- fan/LPC speed constraint  (shaft_lp = Shaft(epsl, Gearf) constructed above Newton loop)
            trf = sqrt(Tt2 / Tref)
            trl = sqrt(Tt2ac / Tref)
            r1, r1_Nf, r1_Nl = shaft_speed_residual(shaft_lp, Nf, Nl, trf, trl)
            res[1, 1] = r1
            a[1, 1] = r1_Nf * Nf_pf
            a[1, 2] = r1_Nl * Nl_pl
            a[1, 3] = 0.0
            a[1, 4] = r1_Nf * Nf_mf
            a[1, 5] = r1_Nl * Nl_ml
            a[1, 6] = 0.0
            a[1, 7] = 0.0
            a[1, 8] = 0.0
            a[1, 9] = 0.0
            rrel[1] = res[1] / Nl

            #-------------------------------------------------------------------------
            #---- fixed corrected mass flow at LPT IGV (vertical-line LPT map)
            res[2, 1], _ = turbine_mb_residual(turb_lp, mblt)
            a[2, 1] = 0.0
            a[2, 2] = mblt_pl
            a[2, 3] = mblt_ph
            a[2, 4] = mblt_mf
            a[2, 5] = mblt_ml
            a[2, 6] = mblt_mh
            a[2, 7] = mblt_Tb
            a[2, 8] = 0.0
            a[2, 9] = mblt_Mi
            rrel[2] = res[2] / mb_lpt_des

            #-------------------------------------------------------------------------
            #---- fixed corrected mass flow at HPT IGV (vertical-line HPT map)
            res[3, 1], _ = turbine_mb_residual(turb_hp, mbht)
            a[3, 1] = 0.0
            a[3, 2] = mbht_pl
            a[3, 3] = mbht_ph
            a[3, 4] = mbht_mf
            a[3, 5] = mbht_ml
            a[3, 6] = mbht_mh
            a[3, 7] = mbht_Tb
            a[3, 8] = 0.0
            a[3, 9] = mbht_Mi
            rrel[3] = res[3, 1] / mb_hpt_des

            #-------------------------------------------------------------------------
            #---- fan nozzle mass flow, choked or unchoked
            mdotf = mf * sqrt(Tref / Tt2) * pt2 / pref
            mdotf_mf = sqrt(Tref / Tt2) * pt2 / pref
            mdotf_pt2 = mf * sqrt(Tref / Tt2) / pref
            res[4, 1], r4_pt7, r4_ht7, r4_st7, r4_Tt7, _, _, _, _, _, _ =
                  nozzle_massflow_residual(nozzle_fan, alpha, nair,
                        pt18, Tt18, ht18, st18, cpt18, Rt18, p0, mdotf)
            a[4, 1] = r4_pt7 * pt7_pf + r4_ht7 * ht7_pf + r4_st7 * st7_pf + r4_Tt7 * Tt7_pf
            a[4, 2] = 0.0
            a[4, 3] = 0.0
            a[4, 4] = mdotf_mf + mdotf_pt2 * pt2_mf +
                      r4_pt7 * pt7_mf + r4_ht7 * ht7_mf + r4_st7 * st7_mf + r4_Tt7 * Tt7_mf
            a[4, 5] = mdotf_pt2 * pt2_ml + r4_pt7 * pt7_ml
            a[4, 6] = 0.0
            a[4, 7] = 0.0
            a[4, 8] = 0.0
            a[4, 9] = mdotf_pt2 * pt2_Mi + r4_pt7 * pt7_Mi
            rrel[4] = res[4, 1] / mdotf

            #-------------------------------------------------------------------------
            #---- core nozzle mass flow, choked or unchoked
            mdotc = (1.0 - fo + ff) * ml * sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotc_ml = (1.0 - fo + ff) * sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotc_fo = -ml * sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotc_ff = ml * sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotc_pt1_9c = (1.0 - fo + ff) * ml * sqrt(Tref / Tt2ac) / pref

            mdotc_pl = mdotc_ff * ff_pl
            mdotc_ph = mdotc_ff * ff_ph
            mdotc_mf = mdotc_ff * ff_mf + mdotc_fo * fo_mf + mdotc_pt1_9c * pt1_9c_mf
            mdotc_ml = mdotc_ff * ff_ml + mdotc_fo * fo_ml + mdotc_pt1_9c * pt1_9c_ml +
                       mdotc_ml
            mdotc_mh = mdotc_ff * ff_mh
            mdotc_Tb = mdotc_ff * ff_Tb
            mdotc_Mi = mdotc_ff * ff_Mi + mdotc_fo * fo_Mi + mdotc_pt1_9c * pt1_9c_Mi

            res[5, 1], r5_pt5, r5_ht5, r5_st5, r5_Tt5,
                  p5_al, T5_al, h5_al, s5_al, cp5_al, R5_al =
                        nozzle_massflow_residual(nozzle_core, lambdap, nair,
                              pt8, Tt8, ht8, st8, cpt8, Rt8, p0, mdotc)
            a[5, 1] = r5_pt5 * pt5_pf + r5_ht5 * ht5_pf + r5_st5 * st5_pf + r5_Tt5 * Tt5_pf
            a[5, 2] = mdotc_pl + r5_pt5 * pt5_pl + r5_ht5 * ht5_pl + r5_st5 * st5_pl + r5_Tt5 * Tt5_pl
            a[5, 3] = mdotc_ph + r5_pt5 * pt5_ph + r5_ht5 * ht5_ph + r5_st5 * st5_ph + r5_Tt5 * Tt5_ph
            a[5, 4] = mdotc_mf + r5_pt5 * pt5_mf + r5_ht5 * ht5_mf + r5_st5 * st5_mf + r5_Tt5 * Tt5_mf
            a[5, 5] = mdotc_ml + r5_pt5 * pt5_ml + r5_ht5 * ht5_ml + r5_st5 * st5_ml + r5_Tt5 * Tt5_ml
            a[5, 6] = mdotc_mh + r5_pt5 * pt5_mh + r5_ht5 * ht5_mh + r5_st5 * st5_mh + r5_Tt5 * Tt5_mh
            a[5, 7] = mdotc_Tb + r5_pt5 * pt5_Tb + r5_ht5 * ht5_Tb + r5_st5 * st5_Tb + r5_Tt5 * Tt5_Tb
            a[5, 8] =            r5_pt5 * pt5_Pc + r5_ht5 * ht5_Pc + r5_st5 * st5_Pc + r5_Tt5 * Tt5_Pc
            a[5, 9] = mdotc_Mi + r5_pt5 * pt5_Mi + r5_ht5 * ht5_Mi + r5_st5 * st5_Mi + r5_Tt5 * Tt5_Mi
            # Species chain rule: lambdap varies with Newton variables
            u5safe = max(u5, T(0.02) * sqrt(Rt8 * Tt8))
            for i = 1:nair
                  rho5_al_i = p5_al[i] / (R5 * T5) -
                              (rho5 / T5) * T5_al[i] -
                              (rho5 / R5) * R5_al[i]
                  u5_al_i   = -h5_al[i] / u5safe
                  dr5_dlambda_i = -A5 * (rho5_al_i * u5 + rho5 * u5_al_i)
                  a[5, 2] += dr5_dlambda_i * lamp_pl[i]
                  a[5, 3] += dr5_dlambda_i * lamp_ph[i]
                  a[5, 4] += dr5_dlambda_i * lamp_mf[i]
                  a[5, 5] += dr5_dlambda_i * lamp_ml[i]
                  a[5, 6] += dr5_dlambda_i * lamp_mh[i]
                  a[5, 7] += dr5_dlambda_i * lamp_Tb[i]
                  a[5, 9] += dr5_dlambda_i * lamp_Mi[i]
            end
            rrel[5] = res[5, 1] / mdotc

            #-------------------------------------------------------------------------
            #---- LPC mass flow = HPC mass flow + mass offtake
            mdotl = ml * sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotl_ml = sqrt(Tref / Tt2ac) * pt2ac / pref
            mdotl_pt1_9c = mdotl / pt2ac

            mdotl_mf = mdotl_pt1_9c * pt1_9c_mf
            mdotl_ml = mdotl_pt1_9c * pt1_9c_ml + mdotl_ml
            mdotl_Mi = mdotl_pt1_9c * pt1_9c_Mi


            mdoth = mh * sqrt(Tref / Tt2_5c) * pt2_5c / pref
            mdoth_mh = sqrt(Tref / Tt2_5c) * pt2_5c / pref
            mdoth_Tt2_5c = -0.5 * mdoth / Tt2_5c
            mdoth_pt2_5c = mdoth / pt2_5c

            mdoth_pl = mdoth_Tt2_5c * Tt2_5c_pl + mdoth_pt2_5c * pt2_5c_pl
            mdoth_mf = mdoth_pt2_5c * pt2_5c_mf
            mdoth_ml = mdoth_Tt2_5c * Tt2_5c_ml + mdoth_pt2_5c * pt2_5c_ml
            mdoth_Mi = mdoth_pt2_5c * pt2_5c_Mi

            res[6, 1] = mdoth - mdotl + mofft
            a[6, 1] = 0.0
            a[6, 2] = mdoth_pl
            a[6, 3] = 0.0
            a[6, 4] = mdoth_mf - mdotl_mf
            a[6, 5] = mdoth_ml - mdotl_ml
            a[6, 6] = mdoth_mh
            a[6, 7] = 0.0
            a[6, 8] = 0.0
            a[6, 9] = mdoth_Mi - mdotl_Mi
            rrel[6] = res[6, 1] / mdoth

            #-------------------------------------------------------------------------

            if opt_calc_call == CalcMode.FixedTt4OffDes
                  #----- specified Tt4 constraint
                  res[7, 1] = Tt4 - Tt4spec
                  a[7, 1] = 0.0
                  a[7, 2] = 0.0
                  a[7, 3] = 0.0
                  a[7, 4] = 0.0
                  a[7, 5] = 0.0
                  a[7, 6] = 0.0
                  a[7, 7] = Tt4_Tb
                  a[7, 8] = 0.0
                  a[7, 9] = 0.0
                  rrel[7] = res[7, 1] / Tt4spec

                  #--------------------------------------------------------
            else
                  #----- specified thrust constraint

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- fan nozzle flow 7-8, use alpha mass fraction (air)
                  pfn = p0 / pt18
                  pfn_pt7 = -pfn / pt18
                  epi = 1.0
                  p8, T8, h8, s8, cp8, R8,
                  p8_pt7,
                  p8_st7, T8_st7, h8_st7, s8_st7,
                  p8_pfn, T8_pfn, h8_pfn, s8_pfn,
                  p8_epi, T8_epi, h8_epi, s8_epi,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_pratd(alpha, nair,
                        pt18, Tt18, ht18, st18, cpt18, Rt18, pfn, epi)
                  h8_pt7 = h8_pfn * pfn_pt7

                  h8_pf = h8_pt7 * pt7_pf + h8_st7 * st7_pf
                  h8_mf = h8_pt7 * pt7_mf + h8_st7 * st7_mf
                  h8_ml = h8_pt7 * pt7_ml
                  h8_Mi = h8_pt7 * pt7_Mi

                  if (ht18 > h8)
                        u8 = sqrt(2.0 * (ht18 - h8))
                        u8_pf = (ht7_pf - h8_pf) / u8
                        u8_mf = (ht7_mf - h8_mf) / u8
                        u8_ml = (ht7_ml - h8_ml) / u8
                        u8_Mi = (ht7_Mi - h8_Mi) / u8
                  else
                        u8 = 0.0
                        u8tmp = max(0.2 * u0, 0.01 * sqrt(Rt18 * Tt18))
                        u8_pf = (ht7_pf - h8_pf) / u8tmp
                        u8_mf = (ht7_mf - h8_mf) / u8tmp
                        u8_ml = (ht7_ml - h8_ml) / u8tmp
                        u8_Mi = (ht7_Mi - h8_Mi) / u8tmp
                  end
                  rho8 = p8 / (R8 * T8)

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- core nozzle flow 5-6, use lambda mass fraction (combustion products)
                  pcn = p0 / pt8
                  pcn_pt5 = -pcn / pt8
                  epi = 1.0
                  p6, T6, h6, s6, cp6, R6,
                  p6_pt5,
                  p6_st5, T6_st5, h6_st5, s6_st5,
                  p6_pcn, T6_pcn, h6_pcn, s6_pcn,
                  p6_epi, T6_epi, h6_epi, s6_epi,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_pratd(lambdap, nair,
                        pt8, Tt8, ht8, st8, cpt8, Rt8, pcn, epi)
                  h6_pt5 = h6_pcn * pcn_pt5

                  h6_pf = h6_pt5 * pt5_pf + h6_st5 * st5_pf
                  h6_pl = h6_pt5 * pt5_pl + h6_st5 * st5_pl
                  h6_ph = h6_pt5 * pt5_ph + h6_st5 * st5_ph
                  h6_mf = h6_pt5 * pt5_mf + h6_st5 * st5_mf
                  h6_ml = h6_pt5 * pt5_ml + h6_st5 * st5_ml
                  h6_mh = h6_pt5 * pt5_mh + h6_st5 * st5_mh
                  h6_Tb = h6_pt5 * pt5_Tb + h6_st5 * st5_Tb
                  h6_Pc = h6_pt5 * pt5_Pc + h6_st5 * st5_Pc
                  h6_Mi = h6_pt5 * pt5_Mi + h6_st5 * st5_Mi
                  for i = 1:nair
                        h6_pl = h6_pl + h_al[i] * lamp_pl[i]
                        h6_ph = h6_ph + h_al[i] * lamp_ph[i]
                        h6_mf = h6_mf + h_al[i] * lamp_mf[i]
                        h6_ml = h6_ml + h_al[i] * lamp_ml[i]
                        h6_mh = h6_mh + h_al[i] * lamp_mh[i]
                        h6_Tb = h6_Tb + h_al[i] * lamp_Tb[i]
                        h6_Mi = h6_Mi + h_al[i] * lamp_Mi[i]
                  end

                  if (ht8 > h6)
                        u6 = sqrt(2.0 * (ht8 - h6))
                        u6_pf = (ht5_pf - h6_pf) / u6
                        u6_pl = (ht5_pl - h6_pl) / u6
                        u6_ph = (ht5_ph - h6_ph) / u6
                        u6_mf = (ht5_mf - h6_mf) / u6
                        u6_ml = (ht5_ml - h6_ml) / u6
                        u6_mh = (ht5_mh - h6_mh) / u6
                        u6_Tb = (ht5_Tb - h6_Tb) / u6
                        u6_Pc = (ht5_Pc - h6_Pc) / u6
                        u6_Mi = (ht5_Mi - h6_Mi) / u6
                  else
                        u6 = 0.0
                        u6tmp = max(0.2 * u0, 0.01 * sqrt(Rt8 * Tt8))
                        u6_pf = (ht5_pf - h6_pf) / u6tmp
                        u6_pl = (ht5_pl - h6_pl) / u6tmp
                        u6_ph = (ht5_ph - h6_ph) / u6tmp
                        u6_mf = (ht5_mf - h6_mf) / u6tmp
                        u6_ml = (ht5_ml - h6_ml) / u6tmp
                        u6_mh = (ht5_mh - h6_mh) / u6tmp
                        u6_Tb = (ht5_Tb - h6_Tb) / u6tmp
                        u6_Mi = (ht5_Mi - h6_Mi) / u6tmp
                  end

                  BPR, BPR_mf_dir, BPR_ml_dir, BPR_pt2, BPR_pt1_9c =
                      bypass_ratio(splitter, mf, ml, pt2, pt2ac, Tt2, Tt2ac)
                  BPR_mf = BPR_mf_dir + BPR_pt2 * pt2_mf + BPR_pt1_9c * pt1_9c_mf
                  BPR_ml = BPR_ml_dir + BPR_pt2 * pt2_ml + BPR_pt1_9c * pt1_9c_ml
                  BPR_Mi =              BPR_pt2 * pt2_Mi  + BPR_pt1_9c * pt1_9c_Mi

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- overall thrust
                  if (u0 == 0.0)
                        Finl = 0.0
                  else
                        Finl = Phiinl / u0
                  end
                  F = ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9) * mdotl + Finl
                  F_ff = u6 * mdotl
                  F_u6 = (1.0 - fo + ff) * mdotl
                  F_BPR = (u8 - u0) * mdotl
                  F_u8 = BPR * mdotl
                  F_fo = (-u6 + u9) * mdotl
                  F_mdotl = (1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9

                  F_pf = F_u6 * u6_pf + F_u8 * u8_pf
                  F_pl = F_ff * ff_pl + F_u6 * u6_pl
                  F_ph = F_ff * ff_ph + F_u6 * u6_ph
                  F_mf = F_ff * ff_mf + F_u6 * u6_mf + F_u8 * u8_mf + F_BPR * BPR_mf
                  F_ml = F_ff * ff_ml + F_u6 * u6_ml + F_u8 * u8_ml + F_BPR * BPR_ml
                  F_mh = F_ff * ff_mh + F_u6 * u6_mh
                  F_Tb = F_ff * ff_Tb + F_u6 * u6_Tb
                  F_Pc = F_u6 * u6_Pc
                  F_Mi = F_ff * ff_Mi + F_u6 * u6_Mi + F_u8 * u8_Mi + F_BPR * BPR_Mi

                  F_mf = F_mf + F_fo * fo_mf + F_mdotl * mdotl_mf
                  F_ml = F_ml + F_fo * fo_ml + F_mdotl * mdotl_ml
                  F_Mi = F_Mi + F_fo * fo_Mi + F_mdotl * mdotl_Mi

                  #----- specified-F constraint
                  res[7, 1] = F - Fspec
                  a[7, 1] = F_pf
                  a[7, 2] = F_pl
                  a[7, 3] = F_ph
                  a[7, 4] = F_mf
                  a[7, 5] = F_ml
                  a[7, 6] = F_mh
                  a[7, 7] = F_Tb
                  a[7, 8] = F_Pc
                  a[7, 9] = F_Mi
                  rrel[7] = res[7, 1] / max(Fspec, F, 1e-6)


                  #      res[7] = Tt8
                  #      a(7,1) = Tt5_pf
                  #      a(7,2) = Tt5_pl
                  #      a(7,3) = Tt5_ph
                  #      a(7,4) = Tt5_mf
                  #      a(7,5) = Tt5_ml
                  #      a(7,6) = Tt5_mh
                  #      a(7,7) = Tt5_Tb
                  #      a(7,8) = Tt5_Pc
                  #      a(7,9) = Tt5_Mi


            end


            #------------------------------------------------------
            #---- LPT power matched to fan+LPC power
            epi = 1.0 / eplt
            epi_eplt = -epi / eplt
            pt5h, Tt5h, ht5h, st5h, cpt5h, Rt5h,
            pt5h_st4_5,
            pt5h_pt4_5,
            pt5h_epi,
            pt5h_ht4_5, Tt5h_ht4_5, ht5h_ht4_5, st5h_ht4_5,
            pt5h_dhlt, Tt5h_dhlt, ht5h_dhlt, st5h_dhlt,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_delhd(lambdap, nair,
                  pt4_5, Tt4_5, ht4_5, st4_5, cpt4_5, Rt4_5, dhlt, epi)

            pt5h_eplt = pt5h_epi * epi_eplt

            pt5h_pf = pt5h_eplt * eplt_pf +
                      pt5h_dhlt * dhlt_pf
            pt5h_pl = pt5h_st4_5 * st4_5_pl +
                      pt5h_pt4_5 * pt4_5_pl + pt5h_eplt * eplt_pl +
                      pt5h_ht4_5 * ht4_5_pl + pt5h_dhlt * dhlt_pl
            pt5h_ph = pt5h_st4_5 * st4_5_ph +
                      pt5h_pt4_5 * pt4_5_ph + pt5h_eplt * eplt_ph +
                      pt5h_ht4_5 * ht4_5_ph + pt5h_dhlt * dhlt_ph
            pt5h_mf = pt5h_st4_5 * st4_5_mf +
                      pt5h_pt4_5 * pt4_5_mf + pt5h_eplt * eplt_mf +
                      pt5h_ht4_5 * ht4_5_mf + pt5h_dhlt * dhlt_mf
            pt5h_ml = pt5h_st4_5 * st4_5_ml +
                      pt5h_pt4_5 * pt4_5_ml + pt5h_eplt * eplt_ml +
                      pt5h_ht4_5 * ht4_5_ml + pt5h_dhlt * dhlt_ml
            pt5h_mh = pt5h_st4_5 * st4_5_mh +
                      pt5h_pt4_5 * pt4_5_mh + pt5h_eplt * eplt_mh +
                      pt5h_ht4_5 * ht4_5_mh + pt5h_dhlt * dhlt_mh
            pt5h_Tb = pt5h_st4_5 * st4_5_Tb +
                      pt5h_pt4_5 * pt4_5_Tb + pt5h_eplt * eplt_Tb +
                      pt5h_ht4_5 * ht4_5_Tb + pt5h_dhlt * dhlt_Tb
            pt5h_Mi = pt5h_st4_5 * st4_5_Mi +
                      pt5h_pt4_5 * pt4_5_Mi + pt5h_eplt * eplt_Mi +
                      pt5h_ht4_5 * ht4_5_Mi + pt5h_dhlt * dhlt_Mi

            for i = 1:nair
                  pt5h_pl = pt5h_pl + p_al[i] * lamp_pl[i]
                  pt5h_ph = pt5h_ph + p_al[i] * lamp_ph[i]
                  pt5h_mf = pt5h_mf + p_al[i] * lamp_mf[i]
                  pt5h_ml = pt5h_ml + p_al[i] * lamp_ml[i]
                  pt5h_mh = pt5h_mh + p_al[i] * lamp_mh[i]
                  pt5h_Tb = pt5h_Tb + p_al[i] * lamp_Tb[i]
                  pt5h_Mi = pt5h_Mi + p_al[i] * lamp_Mi[i]
            end

            res[8, 1] = pt5 - pt5h
            a[8, 1] = pt4_9_pf - pt5h_pf
            a[8, 2] = pt4_9_pl - pt5h_pl
            a[8, 3] = pt4_9_ph - pt5h_ph
            a[8, 4] = pt4_9_mf - pt5h_mf
            a[8, 5] = pt4_9_ml - pt5h_ml
            a[8, 6] = pt4_9_mh - pt5h_mh
            a[8, 7] = pt4_9_Tb - pt5h_Tb
            a[8, 8] = pt4_9_Pc
            a[8, 9] = pt4_9_Mi - pt5h_Mi
            rrel[8] = res[8, 1] / pt5
            # ===========================================================================
            #---- inlet Mach - area constraint
            mfA = mf * sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
            mfA_mf = mfA * (pt2_mf / pt2 - rho2_mf / rho2) +
                     sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
            mfA_ml = mfA * (pt2_ml / pt2 - rho2_ml / rho2)
            mfA_Mi = mfA * (pt2_Mi / pt2 - rho2_Mi / rho2 - u2_Mi / u2)

            mlA = ml * sqrt(Tref / Tt2ac) * pt2ac / pref / (rho1_9c * u1_9c)
            mlA_mf = mlA * (pt1_9c_mf / pt2ac - rho1_9c_mf / rho1_9c)
            mlA_ml = mlA * (pt1_9c_ml / pt2ac - rho1_9c_ml / rho1_9c) +
                     sqrt(Tt2ac / Tref) * pref / pt2ac / (rho1_9c * u1_9c)
            mlA_Mi = mlA * (pt1_9c_Mi / pt2ac - rho1_9c_Mi / rho1_9c - u1_9c_Mi / u1_9c)

            res[9, 1] = mfA + mlA - A2
            a[9, 1] = 0.0
            a[9, 2] = 0.0
            a[9, 3] = 0.0
            a[9, 4] = mfA_mf + mlA_mf
            a[9, 5] = mfA_ml + mlA_ml
            a[9, 6] = 0.0
            a[9, 7] = 0.0
            a[9, 8] = 0.0
            a[9, 9] = mfA_Mi + mlA_Mi
            rrel[9] = res[9, 1]


            # ===========================================================================

            for i = 1:9
                  rsav[i] = res[i, 1]
                  for k = 1:9
                        asav[i, k] = a[i, k]
                  end
            end


            if (iter == -2)
                  for i = 1:9
                        a_dlls[i, 1] = a[i, 1]
                        a_dlls[i, 2] = a[i, 2]
                        a_dlls[i, 3] = a[i, 3]
                        a_dlls[i, 4] = a[i, 4]
                        a_dlls[i, 5] = a[i, 5]
                        a_dlls[i, 6] = a[i, 6]
                        a_dlls[i, 7] = a[i, 7]
                        a_dlls[i, 8] = a[i, 8]
                        a_dlls[i, 9] = a[i, 9]
                        res_dlls[i] = res[i, 1]
                  end
            elseif (iter == -1)

                  for i = 1:9
                        dd = (res[i, 1] - res_dlls[i]) / eps
                        aa = (a[i, j] + a_dlls[i, j]) * 0.5
                        compare(ss, aa, dd)
                        println("Ja", i, j, aa, ss, "Jd", i, j, dd, ss)
                  end

                  error("tfoper.jl error due to iter == -1")

            end


            if (iter < 0)

                  println("TFOPER: Convergence failed.  opt_call_calc=", opt_calc_call)

                  Tt4 = Tb
                  pt8 = Pc
                  BPR = mf / ml * sqrt(Tt2ac / Tt2) * pt2 / pt2ac

                  println("pt12 Tt12 =", pt12, Tt12)
                  println("pt2  Tt2  =", pt2, Tt2)
                  println("pt2_5 Tt2_5 =", pt2_5, Tt2_5)
                  println("pt3  Tt3  =", pt3, Tt3)
                  println("pt4  Tt4  =", pt4, Tt4)
                  println("pt4_1 Tt4_1 =", pt4_1, Tt4_1)
                  println("pt4_5 Tt4_5 =", pt4_5, Tt4_5)
                  println("pt8  Tt8  =", pt8, Tt8)
                  println("p5        =", p5)
                  println("FPR  BPR  =", pf, BPR)

                  Lconv = false
                  return TSFC, Fsp, hfuel, ff,
                  Feng, mcore,
                  pif, pilc, pihc,
                  mbf, mblc, mbhc,
                  Nbf, Nblc, Nbhc,
                  Tt0, ht0, pt0, cpt0, Rt0,
                  Tt12, ht12, pt12, cpt12, Rt12,
                  Tt2a, ht2a, pt2a, cpt2a, Rt2a,
                  Tt2, ht2, pt2, cpt2, Rt2,
                  Tt13, ht13, pt13, cpt13, Rt13,
                  Tt2_5, ht2_5, pt2_5, cpt2_5, Rt2_5,
                  Tt3, ht3, pt3, cpt3, Rt3,
                  Tt4, ht4, pt4, cpt4, Rt4,
                  Tt4_1, ht4_1, pt4_1, cpt4_1, Rt4_1,
                  Tt4_5, ht4_5, pt4_5, cpt4_5, Rt4_5,
                  Tt5, ht5, pt5, cpt5, Rt5,
                  Tt8, ht8, pt8, cpt8, Rt8,
                  Tt18, ht18, pt18, cpt18, Rt18,
                  u0,
                  T2, u2, p2, cp2, R2, M2,
                  T2_5c, u2_5c, p2_5c, cp2_5c, R2_5c, M2_5c,
                  T5, u5, p5, cp5, R5, M5,
                  T6, u6, p6, cp6, R6, M6, A6,
                  T7, u7, p7, cp7, R7, M7,
                  T8, u8, p8, cp8, R8, M8, A8,
                  u9, A9,
                  epf, eplc, ephc, epht, eplt,
                  etaf, etalc, etahc, etaht, etalt,
                  Lconv
            end

            #---- solve Newton system and set Newton deltas
            res = gaussn(9, 9, a, res, 1)

            dMi = -res[9, 1]
            if (Mi >= Mimax && dMi > 0.0)
                  #----- Fan face is approaching choking, possibly due to an iteration transient
                  #-     Artificially limit it to Mimax
                  for i = 1:9
                        res[i, 1] = rsav[i]
                        for k = 1:9
                              a[i, k] = asav[i, k]
                        end
                  end

                  i = 9
                  for k = 1:9
                        a[i, k] = 0.0
                  end
                  res[i, 1] = Mi - Mimax
                  a[i, i] = 1.0

                  res = gaussn(9, 9, a, res, 1)
            end


            dpf = -res[1, 1]
            dpl = -res[2, 1]
            dph = -res[3, 1]
            dmf = -res[4, 1]
            dml = -res[5, 1]
            dmh = -res[6, 1]
            dTb = -res[7, 1]
            dPc = -res[8, 1]
            dMi = -res[9, 1]

            #---- max relative change
            dmax = max(abs(dpf) / pf,
                  abs(dpl) / pl,
                  abs(dph) / ph,
                  abs(dmf) / mf,
                  abs(dml) / ml,
                  abs(dmh) / mh,
                  abs(dTb) / Tb,
                  abs(dPc) / Pc,
                  abs(dMi) / Mi)

            #---- max,min allowed changes 
            pf0 = 1.0
            pl0 = 1.0

            dpfmax = 0.30 * (pf - pf0)
            dplmax = 0.25 * (pl - pl0)
            dphmax = 0.25 * (ph - 1.0)
            dmfmax = 0.20 * mf
            dmlmax = 0.20 * ml
            dmhmax = 0.20 * mh
            dTbmax = 0.50 * (Tb - Tt3)
            dPcmax = 1.00 * (Pc - p0)
            dMimax = 1.001 * (Mimax - Mi)

            dpfmin = -0.30 * (pf - pf0)
            dplmin = -0.25 * (pl - pl0)
            dphmin = -0.25 * (ph - 1.0)
            dmfmin = -0.20 * mf
            dmlmin = -0.20 * ml
            dmhmin = -0.20 * mh
            dTbmin = -0.30 * (Tb - Tt3)
            dPcmin = -0.50 * (Pc - p0)
            dMimin = -0.20 * Mi

            #---- set underrelaxation factor, 
            #-    if needed to limit any one change to its max value
            rlx = 1.0
            
            if (rlx * dpf > dpfmax)
                  rlx = dpfmax / dpf
            end
            if (rlx * dpl > dplmax)
                  rlx = dplmax / dpl
            end
            if (rlx * dph > dphmax)
                  rlx = dphmax / dph
            end
            if (rlx * dmf > dmfmax)
                  rlx = dmfmax / dmf
            end
            if (rlx * dml > dmlmax)
                  rlx = dmlmax / dml
            end
            if (rlx * dmh > dmhmax)
                  rlx = dmhmax / dmh
            end
            if (rlx * dTb > dTbmax)
                  rlx = dTbmax / dTb
            end
            if (rlx * dPc > dPcmax)
                  rlx = dPcmax / dPc
            end
            if (rlx * dMi > dMimax)
                  rlx = dMimax / dMi
            end

            if (rlx * dpf < dpfmin)
                  rlx = dpfmin / dpf
            end
            if (rlx * dpl < dplmin)
                  rlx = dplmin / dpl
            end
            if (rlx * dph < dphmin)
                  rlx = dphmin / dph
            end
            if (rlx * dmf < dmfmin)
                  rlx = dmfmin / dmf
            end
            if (rlx * dml < dmlmin)
                  rlx = dmlmin / dml
            end
            if (rlx * dmh < dmhmin)
                  rlx = dmhmin / dmh
            end
            if (rlx * dTb < dTbmin)
                  rlx = dTbmin / dTb
            end
            if (rlx * dPc < dPcmin)
                  rlx = dPcmin / dPc
            end
            if (rlx * dMi < dMimin)
                  rlx = dMimin / dMi
            end



            rlx = 1.0
            vrlx = " "

            if (rlx * dpf > dpfmax)
                  rlx = dpfmax / dpf
                  vrlx = "pf"
            end
            if (rlx * dpl > dplmax)
                  rlx = dplmax / dpl
                  vrlx = "pl"
            end
            if (rlx * dph > dphmax)
                  rlx = dphmax / dph
                  vrlx = "ph"
            end
            if (rlx * dmf > dmfmax)
                  rlx = dmfmax / dmf
                  vrlx = "mf"
            end
            if (rlx * dml > dmlmax)
                  rlx = dmlmax / dml
                  vrlx = "ml"
            end
            if (rlx * dmh > dmhmax)
                  rlx = dmhmax / dmh
                  vrlx = "mh"
            end
            if (rlx * dTb > dTbmax)
                  rlx = dTbmax / dTb
                  vrlx = "Tb"
            end
            if (rlx * dPc > dPcmax)
                  rlx = dPcmax / dPc
                  vrlx = "Pc"
            end
            if (rlx * dMi > dMimax)
                  rlx = dMimax / dMi
                  vrlx = "Mi"
            end

            if (rlx * dpf < dpfmin)
                  rlx = dpfmin / dpf
                  vrlx = "pf"
            end
            if (rlx * dpl < dplmin)
                  rlx = dplmin / dpl
                  vrlx = "pl"
            end
            if (rlx * dph < dphmin)
                  rlx = dphmin / dph
                  vrlx = "ph"
            end
            if (rlx * dmf < dmfmin)
                  rlx = dmfmin / dmf
                  vrlx = "mf"
            end
            if (rlx * dml < dmlmin)
                  rlx = dmlmin / dml
                  vrlx = "ml"
            end
            if (rlx * dmh < dmhmin)
                  rlx = dmhmin / dmh
                  vrlx = "mh"
            end
            if (rlx * dTb < dTbmin)
                  rlx = dTbmin / dTb
                  vrlx = "Tb"
            end
            if (rlx * dPc < dPcmin)
                  rlx = dPcmin / dPc
                  vrlx = "Pc"
            end
            if (rlx * dMi < dMimin)
                  rlx = dMimin / dMi
                  vrlx = "Mi"
            end

            #If iter>10, a limit cycle may have been reached
            #Apply a relaxation factor that does not oscillate
            rlx_it = 1.0 - 0.6*rand() * (iter > 10) #Relaxation based on RNG
            rlx = rlx * rlx_it

            #---- exit if convergence test is met or if max iterations reached
            if (dmax < toler) | (iter == itmax)

                  # ===============================================================
                  #---- pick up here if converged normally

                  #---- fan nozzle flow 7-8, use alpha mass fraction (air)
                  pt19 = pt18
                  ht19 = ht18
                  Tt19 = Tt18
                  st19 = st18
                  cpt19 = cpt18
                  Rt19 = Rt18

                  pfn = p0 / pt19
                  p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair,
                        pt19, Tt19, ht19, st19, cpt19, Rt19, pfn, 1.0)
                  if (ht19 > h8)
                        u8 = sqrt(2.0 * (ht19 - h8))
                  else
                        u8 = 0.0
                  end
                  rho8 = p8 / (R8 * T8)
                  #
                  #---------------------------------------------------------------
                  #---- core nozzle flow 5-6, use lambda mass fraction (combustion products)
                  pt9 = pt8
                  ht9 = ht8
                  Tt9 = Tt8
                  st9 = st8
                  cpt9 = cpt8
                  Rt9 = Rt8

                  pcn = p0 / pt9
                  p6, T6, h6, s6, cp6, R6 = gas_prat(lambdap, nair,
                        pt9, Tt9, ht9, st9, cpt9, Rt9, pcn, 1.0)
                  if (ht9 > h6)
                        u6 = sqrt(2.0 * (ht9 - h6))
                  else
                        u6 = 0.0
                  end
                  rho6 = p6 / (R6 * T6)

                  # ===============================================================
                  #---- set final corrected mass flows, pressure ratios, corrected speeds
                  mbf = mf
                  mblc = ml
                  mbhc = mh

                  pif = pf
                  pilc = pl
                  pihc = ph

                  Nbf = Nl / Gearf * sqrt(Tt2ac / Tt2)
                  Nblc = Nl
                  Nbhc = Nh

                  Tt4 = Tb
                  pt8 = Pc

                  #---- final core mass flow and bypass ratio
                  mcore = ml * sqrt(Tref / Tt2ac) * pt2ac / pref
                  BPR = mf / ml * sqrt(Tt2ac / Tt2) * pt2 / pt2ac

                  #---- offtake mass ratio
                  fo = mofft / mcore

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                  #---- effective fuel heating value, from states 3, 4  (just for info)
                  cpa = 0.5 * (cpt3 + cpt4)
                  hfuel = cpa * (Tt4 - Tt3 + ffb * (Tt4 - Ttf)) / (etab * ffb)

                  #---- effective fuel heating value, from states 3, 4.1  (just for info)
                  #      cpa = 0.5*(cpt3+cpt4_1)
                  #      hfuel = cpa*((1.0-fo)*(Tt4_1-Tt3) + ff*(Tt4_1-Ttf)) / (etab*ff)
                  #- - - - - - - - - - - - - - - - - - - - - - - - - - -

                  #---- overall effective thrust, including inlet defect credit
                  if (u0 == 0.0)
                        Finl = 0.0
                  else
                        Finl = Phiinl / u0
                  end
                  Feng = ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9) * mcore + Finl

                  #---- overall Fsp and TSFC
                  if (u0 == 0.0)
                        Fsp = 0.0
                  else
                        Fsp = Feng / (u0 * mcore * (1.0 + BPR))
                  end

                  if (Feng == 0.0)
                        TSFC = 0.0
                  else
                        TSFC = (gee * ff * mcore) / Feng
                  end

                  #-------------------------------------------------------------
                  #---- plume areas
                  A8 = BPR * mcore / (rho8 * u8)
                  A6 = (1.0 - fo + ff) * mcore / (rho6 * u6)

                  if (u9 == 0.0)
                        A9 = 0.0
                  else
                        A9 = fo * mcore / (rho9 * u9)
                  end

                  #---- plume Mach numbers
                  M8 = u8 / sqrt(T8 * R8 * cp8 / (cp8 - R8))
                  M6 = u6 / sqrt(T6 * R6 * cp6 / (cp6 - R6))

                  #--------------------------------------------------------------
                  #---- calculate component efficiencies  (informative only -- not needed here)
                  etaf = 0.0
                  etalc = 0.0
                  etahc = 0.0
                  etaht = 0.0
                  etalt = 0.0

                  #---- fan
                  pt2_1i, Tt2_1i, ht2_1i, st2_1i, cpt2_1i, Rt2_1i = gas_prat(alpha, nair,
                        pt2, Tt2, ht2, st2, cpt2, Rt2, pif, 1.0)
                  etaf = (ht2_1i - ht2) / (ht13 - ht2)

                  #---- LP compressor
                  pt2_5i, Tt2_5i, ht2_5i, st2_5i, cpt2_5i, Rt2_5i = gas_prat(alpha, nair,
                        pt2ac, Tt2ac, ht2ac, st2ac, cpt2ac, Rt2ac, pilc, 1.0)
                  etalc = (ht2_5i - ht2ac) / (ht2_5 - ht2ac)

                  #---- HP compressor

                  pt3i, Tt3i, ht3i, st3i, cpt3i, Rt3i = gas_prat(alpha, nair,
                        pt2_5c, Tt2_5c, ht2_5c, st2_5c, cpt2_5c, Rt2_5c, pihc, 1.0)
                  etahc = (ht3i - ht2_5c) / (ht3 - ht2_5c)

                  #---- HP turbine
                  piht = pt4_5 / pt4_1
                  pt4_5i, Tt4_5i, ht4_5i, st4_5i, cpt4_5i, Rt4_5i = gas_prat(lambdap, nair,
                        pt4_1, Tt4_1, ht4_1, st4_1, cpt4_1, Rt4_1, piht, 1.0)
                  etaht = (ht4_5 - ht4_1) / (ht4_5i - ht4_1)

                  #---- LP turbine
                  pilt = pt5 / pt4_5
                  pt4_9i, Tt4_9i, ht4_9i, st4_9i, cpt4_9i, Rt4_9i = gas_prat(lambdap, nair,
                        pt4_5, Tt4_5, ht4_5, st4_5, cpt4_5, Rt4_5, pilt, 1.0)
                  etalt = (ht8 - ht4_5) / (ht4_9i - ht4_5)

                  #--------------------------------------------------------------
                  #---- calculate fan-face static quantities from Mach number variable Mi
                  M2 = Mi
                  p2, T2, h2, s2, cp2, R2 = gas_mach(alpha, nair,
                        pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, M2, 1.0)
                  u2 = sqrt(2.0 * (ht2 - h2))

                  #---- calculate HPC-face static quantities from core mass flow
                  #-    (first check for choking)
                  M2_5s = 0.99
                  p2_5s, T2_5s, h2_5s, s2_5s, cp2_5s, R2_5s = gas_mach(alpha, nair,
                        pt2_5c, Tt2_5c, ht2_5c, st2_5c, cpt2_5c, Rt2_5c, 0.0, M2_5s, 1.0)
                  u2_5s = sqrt(2.0 * (ht2_5c - h2_5s))
                  mdot2_5s = p2_5s / (R2_5s * T2_5s) * u2_5s * A2_5

                  if ((1.0 - fo) * mcore >= mdot2_5s)
                        #----- HPC face is choked... artificially use sonic quantities
                        u2_5c = u2_5s
                        p2_5c = p2_5s
                        T2_5c = T2_5s
                        cp2_5c = cp2_5s
                        R2_5c = R2_5s
                        M2_5c = M2_5s
                  else
                        #----- normal case... set static quantities from prescribed mass flux
                        Mguess = min(M2_5, 0.90)
                        mp2_5 = (1.0 - fo) * mcore / A2_5
                        p2_5c, T2_5c, h2_5c, s2_5c, cp2_5c, R2_5c = gas_mass(alpha, nair,
                              pt2_5c, Tt2_5c, ht2_5c, st2_5c, cpt2_5c, Rt2_5c, mp2_5, Mguess)
                        u2_5c = sqrt(2.0 * (ht2_5c - h2_5c))
                        M2_5c = u2_5c / sqrt(T2_5c * R2_5c * cp2_5c / (cp2_5c - R2_5c))
                  end

                  if (ht8 < h5) #error if total enthalpy at nozzle inlet is lower than static enthalpy
                        error("ht8 < h5 : ", ht8, h5)
                  end
                  
                  if (dmax < toler) #Convergence achieved
                        Lconv = true

                  else #Finish run as max iterations reached
                        Lconv = false
                        
                        # println("TFOPER: Convergence failed.  iTFspec=", iTFspec)

                        # Tt4 = Tb
                        # pt8 = Pc
                        # BPR = mf / ml * sqrt(Tt2ac / Tt2) * pt2 / pt2ac

                        # println("pt12 Tt12 =", pt12, Tt12)
                        # println("pt2  Tt2  =", pt2, Tt2)
                        # println("pt2_5 Tt2_5 =", pt2_5, Tt2_5)
                        # println("pt3  Tt3  =", pt3, Tt3)
                        # println("pt4  Tt4  =", pt4, Tt4)
                        # println("pt4_1 Tt4_1 =", pt4_1, Tt4_1)
                        # println("pt4_5 Tt4_5 =", pt4_5, Tt4_5)
                        # println("pt8  Tt8  =", pt8, Tt8)
                        # println("p5        =", p5)
                        # println("FPR  BPR  =", pf, BPR)
                  end
                 
                  return TSFC, Fsp, hfuel, ff,
                  Feng, mcore,
                  pif, pilc, pihc,
                  mbf, mblc, mbhc,
                  Nbf, Nblc, Nbhc,
                  Tt0, ht0, pt0, cpt0, Rt0,
                  Tt12, ht12, pt12, cpt12, Rt12,
                  Tt2a, ht2a, pt2a, cpt2a, Rt2a,
                  Tt2ac, ht2ac, pt2ac, cpt2ac, Rt2ac,
                  Tt2, ht2, pt2, cpt2, Rt2,
                  Tt13, ht13, pt13, cpt13, Rt13,
                  Tt2_5, ht2_5, pt2_5, cpt2_5, Rt2_5,
                  Tt2_5c, ht2_5c, pt2_5c, cpt2_5c, Rt2_5c,
                  Tt3, ht3, pt3, cpt3, Rt3,
                  Tt4, ht4, pt4, cpt4, Rt4,
                  Tt4_1, ht4_1, pt4_1, cpt4_1, Rt4_1,
                  Tt4_5, ht4_5, pt4_5, cpt4_5, Rt4_5,
                  Tt5, ht5, pt5, cpt5, Rt5,
                  Tt8, ht8, pt8, cpt8, Rt8,
                  Tt18, ht18, pt18, cpt18, Rt18,
                  u0,
                  T2, u2, p2, cp2, R2, M2,
                  T2_5c, u2_5c, p2_5c, cp2_5c, R2_5c, M2_5c,
                  T5, u5, p5, cp5, R5, M5,
                  T6, u6, p6, cp6, R6, M6, A6,
                  T7, u7, p7, cp7, R7, M7,
                  T8, u8, p8, cp8, R8, M8, A8,
                  u9, A9,
                  epf, eplc, ephc, epht, eplt,
                  etaf, etalc, etahc, etaht, etalt,
                  Lconv

            end
            

            #---- Newton update
            pf = pf + rlx * dpf
            pl = pl + rlx * dpl
            ph = ph + rlx * dph
            mf = mf + rlx * dmf
            ml = ml + rlx * dml
            mh = mh + rlx * dmh
            Tb = Tb + rlx * dTb
            Pc = Pc + rlx * dPc
            Mi = Mi + rlx * dMi

            Mi = min(Mi, Mimax)

      end # with next Newton iteration

end # tfoper


