"""
      ductedfansize!(gee, M0, T0, p0, a0, M2,
            Feng, Phiinl, Kinl, iBLIc,
            pif,
            pid, pifn, 
            epf0
            )

Ducted fan performance and sizing routine.
      
This model is based on the turbofan model in `tfsize()`, stripped of the core.
      
The gas routines reside in the following source files:
    gascalc.f  Routines for various processes (compressor, turbine, combustor, etc)
    gasfun.f   Routines for computing cp[T], h[t], sigma[T], R, called by the routines in gascalc.f
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `M2`:      fan-face Mach number
    - `Feng`:    required net thrust  (PK_inl+PK_out-Phi_jet)/u0  =  sum( mdot u)
    - `Phiinl`:  inlet ingested dissipation
    - `pif`:     fan      pressure ratio  ( = pt7 /pt2)
    - `pid`:     diffuser pressure ratio  ( = pt2 /pt0)
    - `pib`:     burner   pressure ratio  ( = pt4 /pt3)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    fan max polytropic efficiency

      **Outputs:**
    - `Fsp`:     specific thrust  = F / (mdot u0) = F / ((1+BPR) mdot_core u0)
    - `Pfan`:    power required to drive fan [W]
    - `mcore`:   core mass flow = mdot_core  [kg/s]
    - `A2`:      fan-face area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `A8`:      fan  plume  area [m^2]
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= dh/dT)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:      static temperature
    - `u?`:      velocity
    - `epf`:     fan polytropic efficiency
    - `etaf`:    fan overall efficiency
    - `Lconv`:   T if convergence was successful, F otherwise

    The "?" symbol denotes the station index:
      0  freestream
      18 fan face outside of casing BLs
      19 fan face over LPC portion
      2  fan face over fan portion
      21 fan exit
      7  fan nozzle
      8  fan flow downstream
"""
function ductedfansize!(gee, M0, T0, p0, a0, M2,
      Feng, Phiinl, Kinl, iBLIc,
      pif,
      pid, pifn, 
      epf0,
      Δh_radiator,
      Δp_radiator
      )

      n = 6

      # from 'airfrac.inc'
      # air fractions  
      #        N2      O2      CO2    H2O      Ar       fuel
      alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]

      #---- fractional core mass flow convergence tolerance
      toler = 1.0e-12

      #---- number of air constitutents (all but fuel)
      nair = n - 1

      mfan = 0.0 #Initialize

      # ===============================================================
      #---- freestream static quantities
      s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0)
      gam0 = cp0 / (cp0 - R0)
      u0 = M0 * a0

      # ===============================================================
      #---- freestream total quantities
      hspec = h0 + 0.5 * u0^2
      Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
      Tt0 = gas_tset(alpha, nair, hspec, Tguess)

      st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
      pt0 = p0 * exp((st0 - s0) / Rt0)
      at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))

      # ===============================================================
      #---- diffuser flow 0-2 via shared Inlet component (BLI parameters included)
      _inl  = Inlet(pid; Kinl=Float64(Kinl), eng_has_BLI_cores=(iBLIc != 0))
      _fs0  = FlowStation{Float64}(Tt0, ht0, pt0, cpt0, Rt0,
                  SVector{5,Float64}(alpha[1], alpha[2], alpha[3], alpha[4], alpha[5]))
      _fs0.st = st0  # entropy-complement not accepted by the 6-arg constructor
      _fs18 = FlowStation{Float64}()
      inlet_diffuser!(_fs18, _fs0, _inl)
      Tt18, ht18, pt18, cpt18, Rt18, st18 =
          _fs18.Tt, _fs18.ht, _fs18.pt, _fs18.cpt, _fs18.Rt, _fs18.st

      # BLI output stations (ducted fan has no core; _fs19 satisfies interface)
      _fs2  = FlowStation{Float64}()
      _fs19 = FlowStation{Float64}()

      #---- initial guesses for station 2 and 1.9
      pt2 = pt18
      Tt2 = Tt18

      # Fan compressor component — design anchors: pi_fan_des = pif (sizing point IS design point),
      # mb_fan_des = 1.0 (corrected flow normalised to design), NbD = 1.0.
      _comp_fan = Compressor(pif, 1.0, 1.0, epf0, 0.60, FanMap)

      npass = 60

      for ipass = 1:npass

            # ===============================================================
            #---- set fan inlet conditions corrected for BLI
            if (ipass == 1)
                  #----- mfan not yet computed; initialise st2 = st18 (sbfan = 0)
                  _fs2.Tt    = _fs18.Tt;  _fs2.ht    = _fs18.ht;  _fs2.pt   = _fs18.pt
                  _fs2.cpt   = _fs18.cpt; _fs2.Rt    = _fs18.Rt;  _fs2.st   = _fs18.st
                  _fs2.alpha = _fs18.alpha
            else
                  #----- account for inlet BLI defect via mass-averaged entropy
                  #      corrected-flow normalisation matches ductedfanoper!/inlet_bli_mixing!
                  mf_corr = mfan * sqrt(Tt2 / Tref) / (pt2 / pref)
                  inlet_bli_mixing!(_fs2, _fs19, _fs18, _fs0, _inl,
                                    mf_corr, 0.0, M2, at0, gam0, Tref, pref)
            end

            #---- note: BL is assumed adiabatic,
            #-     so Tt2,ht2,st2,cpt2,Rt2  will not change due to BL ingestion
            Tt2  = _fs2.Tt
            ht2  = _fs2.ht
            st2  = _fs2.st
            cpt2 = _fs2.cpt
            Rt2  = _fs2.Rt
            pt2  = _fs2.pt


            # ===============================================================
            #---- fan flow 2-7 via shared Compressor component (efficiency map only;
            #     gas_prat preserves exact FP evaluation order of the original sizing code)
            _, epf, _, _, _, _ = compressor_efficiency(_comp_fan, pif, 1.0)

            pt21, Tt21, ht21, st21, cpt21, Rt21 = gas_prat(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, pif, epf)

            # ===============================================================
            #---- Radiator heat exchanger
            pt7 = pt21 * pifn - Δp_radiator
            ht7 = ht21 + Δh_radiator
      
            Tt7 = gas_tset(alpha, nair, ht7, Tt21)
            st7, _, ht7, _, cpt7, Rt7 = gassum(alpha, nair, Tt7)

            # ===============================================================
            #---- Fan nozzle (7) and plume (8) via shared Nozzle component.
            # pt7 already incorporates pifn and radiator Δp; no additional
            # pressure loss is applied inside nozzle_exit (pn = 1).
            fan_nozzle = Nozzle(1.0, 0.0)
            p7, T7, h7, s7, cp7, R7, u7, rho7, M7 = nozzle_exit(
                  fan_nozzle, alpha, nair, pt7, Tt7, ht7, st7, cpt7, Rt7, p0)

            if (M7 < 1.0)
                  #----- subsonic nozzle: plume fully expanded, same state as nozzle
                  p8, T8, h8, s8, cp8, R8, u8, rho8 = p7, T7, h7, s7, cp7, R7, u7, rho7
            else
                  #----- choked nozzle: expand from pt7 to ambient for fan plume (8)
                  p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair,
                        pt7, Tt7, ht7, st7, cpt7, Rt7, p0/pt7, 1.0)
                  if (h8 >= ht7)
                        error("ductedfansize: negative fan plume velocity",
                              "\n\tpt7, Tt7 = ", pt7, " Pa, ", Tt7, " K",
                              "\n\tp8,  T8  = ", p8,  " Pa, ", T8,  " K")
                  end
                  u8   = sqrt(2.0 * (ht7 - h8))
                  rho8 = p8 / (R8 * T8)
            end

            # ===============================================================
            #---- size mass flow
            #---- store current values for better update, convergence checks
            mfold = mfan

            #---- added effective net thrust from dissipation in ingested streamtube
            if (u0 == 0.0)
                Finl = 0.0
            else
                Finl = Phiinl / u0
            end
            #---- set core mass flow from specified effective net thrust
            mfan = (Feng - Finl) /
                  (u7 + (p7 - p0)/(rho7 * u7) - u0)

            #---- overall Fsp and TSFC
            Fsp = Feng / (u0 * mfan)
            #---- Fan power 
            Pfan = mfan * (ht21 - ht2)

            TSEC = Pfan/Feng

            #---- size fan  nozzle and plume areas
            A7 = mfan / (rho7 * u7)
            A8 = mfan / (rho8 * u8)

            # ===============================================================
            #---- size fan areas
            p2, T2, h2, s2, cp2, R2 = gas_mach(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, M2, 1.0)
            u2 = sqrt(2.0 * (ht2 - h2))
            rho2 = p2 / (R2 * T2)

            A2 = mfan / (rho2 * u2)

            if (ipass >= 2)
                  dmfrac = 1.0 - mfold / mfan

                  if (abs(dmfrac) < toler)

                        # ===============================================================
                        #---- calculate component efficiencies  (informative only -- not needed here)
                        etaf = 0.0

                        #---- fan
                        pt21i, Tt21i, ht21i, st21i, cpt21i, Rt21i = gas_prat(alpha, nair,
                              pt2, Tt2, ht2, st2, cpt2, Rt2, pif, 1.0)
                        etaf = (ht21i - ht2) / (ht21 - ht2)
                        
                        Lconv = true
                        return TSEC, Fsp, Pfan, mfan,
                        Tt0, ht0, pt0, cpt0, Rt0,
                        Tt18, ht18, pt18, cpt18, Rt18,
                        Tt2, ht2, pt2, cpt2, Rt2,
                        Tt21, ht21, pt21, cpt21, Rt21,
                        Tt7, ht7, pt7, cpt7, Rt7,
                        u0,
                        T2, u2, p2, cp2, R2, A2,
                        T7, u7, p7, cp7, R7, A7,
                        T8, u8, p8, cp8, R8, A8,
                        epf,
                        etaf,
                        Lconv
                  end

            end

      end

      if (npass > 1)
            println("DUCTEDFANSIZE: Convergence failed.  dm/m = ", dmfrac)
      end

end # tfsize