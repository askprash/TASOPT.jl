"""
    tfcalc!(wing, engine, parg, para, eng_hx, ip, ifuel, opt_calc_call, opt_cooling, initializes_engine)

Calls on-design sizing function [`tfsize!`](@ref) or off-design analysis function
[`tfoper!`](@ref) for one operating point `ip`.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `eng_hx::EngineState`: per-point typed engine state supplying HX delta inputs
      (hvapcombustor, PreCDeltah/p, InterCDeltah/p, RegenDeltah/p, TurbCDeltah, HXrecircP).
      Populated by `resetHXs`/`HXOffDesign!` before this call. These fields override the
      (hvapcombustor, PreCDeltah/p, etc.) so that HX delta reads come from typed state
      (tasopt-dti / tasopt-w82).
    - `opt_calc_call::CalcMode.T`:
      - `CalcMode.Sizing`: call on-design sizing routine `tfsize!`
      - `CalcMode.FixedTt4OffDes`: call off-design analysis routine `tfoper!` with specified Tt4
      - `CalcMode.FixedFeOffDes`: call off-design analysis routine `tfoper!` with specified net thrust

    - `opt_cooling::CoolingOpt.T`: turbine cooling model
      - `CoolingOpt.NoCooling`: no cooling mass flow
      - `CoolingOpt.FixedCoolingFlowRatio`: cooling flow ratios `epsrow` are inputs; compute `Tmrow`
      - `CoolingOpt.FixedTmetal`: metal temperatures `Tmrow` are inputs; compute `epsrow`
    - `initializes_engine`:
      - `true`: initialize variables for iteration in `tfoper!`
      - `false`: use current variables as initial guesses in `tfoper!`
"""
function tfcalc!(wing, engine, parg::Vector{Float64}, para, eng_hx::EngineState, ip::Int64, ifuel::Int64,
        opt_calc_call::CalcMode.T, opt_cooling::CoolingOpt.T, initializes_engine::Bool)

        # ── ENTRY: alias per-point EngineState (tasopt-j9l.45.16) ────────────────
        # Typed EngineState is the single source of truth; tfcalc! reads all
        # inputs from eng and writes all outputs directly back to eng.
        # HX delta fields are set by resetHXs/HXOffDesign! before this call.
        # Non-pare outputs (st2ac/st2_5c) accumulate across calls.
        eng = eng_hx

        Lprint = false

        if (Lprint)
                println("entering TFCALC", opt_calc_call, opt_cooling, initializes_engine)
        end

        eng_has_BLI_cores = engine.model.has_BLI_cores

        Gearf = parg[igGearf]
        Tmetal = parg[igTmetal]
        neng = parg[igneng]

        mofWpay = parg[igmofWpay]
        mofWMTO = parg[igmofWMTO]
        PofWpay = parg[igPofWpay]
        PofWMTO = parg[igPofWMTO]

        Wpay = parg[igWpay]
        WMTO = parg[igWMTO]

        TSFC   = eng.TSFC
        Fsp    = eng.Fsp
        hfuel  = eng.hfuel
        Tfuel  = eng.Tfuel
        Tt4    = eng.Tt4
        BPR    = eng.BPR
        pif    = eng.pif
        pilc   = eng.pilc
        pihc   = eng.pihc
        pid    = eng.design.pid
        pib    = eng.design.pib
        pifn   = eng.design.pifn
        pitn   = eng.design.pitn
        epolf  = eng.design.epolf
        epollc = eng.design.epollc
        epolhc = eng.design.epolhc
        epolht = eng.design.epolht
        epollt = eng.design.epollt
        etab   = eng.design.etab
        M2     = eng.design.M2
        M2_5    = eng.design.M25
        # Ambient / freestream scalars — read from typed state.
        # Tt0/ht0/pt0/cpt0/Rt0 are computed outputs of tfsize!/tfoper! and are
        # not pre-read here; they are assigned from the returned tuple below.
        M0    = eng.M0
        T0    = eng.T0
        p0    = eng.p0
        a0    = eng.a0
        u0    = eng.st0.u
        rho0  = eng.rho0
        mu0   = eng.mu0
        mcore = eng.st2.mdot
        dTstrk = eng.design.dTstrk
        StA    = eng.design.StA
        Mtexit = eng.design.Mtexit
        M4a = eng.design.M4a
        ruc = eng.design.ruc
        efilm  = eng.design.efilm
        tfilm  = eng.design.tfilm
        epsl   = eng.design.epsl
        epsh   = eng.design.epsh
        epsrow = zeros(ncrowx)
        Tmrow  = zeros(ncrowx)
        hvap   = eng.hvapcombustor

        #Effect of cooling on HPT efficiency
        epht_fc = eng.design.dehtdfc
        fc0     = eng.design.fc0

        #Heat exchanger variables
        Δh_PreC = eng.PreCDeltah
        Δh_InterC = eng.InterCDeltah
        Δh_Regen = eng.RegenDeltah
        Δh_TurbC = eng.TurbCDeltah
        Δp_PreC = eng.PreCDeltap
        Δp_InterC = eng.InterCDeltap
        Δp_Regen = eng.RegenDeltap

        if opt_cooling == CoolingOpt.FixedCoolingFlowRatio
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        epsrow[icrow] = eng.design.epsrow[icrow]
                end
        elseif opt_cooling == CoolingOpt.FixedTmetal
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        #cc      Tmrow[icrow]  = parg[igTmetal]  (Fortran bare-pare equivalent removed)
                        Tmrow[icrow] = parg[igTmetal]
                end
        end

        #- - - - - - - - - - - - - - - - - - - - - - - 
        #---- set BL ingestion parameters
        if (M0 == 0.0)
                #----- no ingestion for static case
                Phiinl = 0.0
                Kinl = 0.0

        else
                #----- assume engine is at TE of fuselage
                DAfsurf = para[iaDAfsurf]
                KAfTE = para[iaKAfTE]
                fBLIf = parg[igfBLIf]

                #----- assume 85% of wing dissipation is on surface
                fDwake = 0.15
                CDAwing = para[iaCDwing] * wing.layout.S
                DAwsurf = CDAwing * (1.0 - fDwake)
                KAwTE = DAwsurf

                fBLIw = parg[igfBLIw]

                #----- set ingested  PKinl-PVinl = Phiinl  for one engine
                Phiinl = 0.5 * rho0 * u0^3 * (DAfsurf * fBLIf + DAwsurf * fBLIw) / neng
                Kinl = 0.5 * rho0 * u0^3 * (KAfTE * fBLIf + KAwTE * fBLIw) / neng
        end
        #- - - - - - - - - - - - - - - - - - - - - - - 

        #---- mass and power offtakes
        Pofft_HX = eng.HXrecircP #power offtakes to drive heat exchanger recirculation per engine
        mofft = (mofWpay * Wpay + mofWMTO * WMTO) / neng
        Pofft = (PofWpay * Wpay + PofWMTO * WMTO) / neng + Pofft_HX

        Tt25off = eng.st25off.Tt
        pt25off = eng.st25off.pt

        #--------------------------------------------------------------------------
        #Engine model convergence
        eng.ConvFail = 0.0

        # #--------------------------------------------------------------------------
        if opt_calc_call == CalcMode.Sizing
                #----- engine sizing case

                # ── ENTRY ────────────────────────────────────────────────────
                # Alias the top-level EngineState built before the if/else.
                # After tfsize! returns we populate every station and the
                # DesignState from local variables, then project back to pare.
                eng_design = eng

                Fe = eng.Fe

                if (Lprint)
                        println("TFSIZE  M0 p0 =", M0, p0)
                        println("     pif pilc =", pif, pilc)
                        println("        Tt4 F =", Tt4, Fe)
                        println("ncrow =", ncrow)
                        println("        alt_m =", para[iaalt])
                end

                result = tfsize!(gee, M0, T0, p0, a0, M2, M2_5,
                        Fe, Phiinl, Kinl, eng_has_BLI_cores,
                        BPR, pif, pilc, pihc,
                        pid, pib, pifn, pitn,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt25off, pt25off, Tt4,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                        Δp_PreC, Δp_InterC, Δp_Regen)

                # Station Mach numbers — needed for ichoke5/ichoke7 at function exit.
                # Computed here from SizingResult static states; not stored in typed state.
                M5 = result.u5 / sqrt(result.T5 * result.R5 * result.cp5 / (result.cp5 - result.R5))
                M6 = result.u6 / sqrt(result.T6 * result.R6 * result.cp6 / (result.cp6 - result.R6))
                M7 = result.u7 / sqrt(result.T7 * result.R7 * result.cp7 / (result.cp7 - result.R7))
                M8 = result.u8 / sqrt(result.T8 * result.R8 * result.cp8 / (result.cp8 - result.R8))

                if (Lprint)
                        println("exited TFSIZE", mofft / result.mcore)
                end

                # ── EXIT: write SizingResult back to EngineState ─────────────
                _update_engine_sizing!(eng_design, result,
                        BPR, Fe, mofft, Pofft, Phiinl, Kinl, neng,
                        pif, pilc, pihc, ruc, M4a, epolht, epollt, Gearf, parg)

                #--------------------------------------------------------------------------
        else
                #----- off-design operation case

                # ── ENTRY ────────────────────────────────────────────────────
                # Alias the top-level EngineState built before the if/else.
                # After tfoper! returns we populate every station from local
                # variables directly into typed EngineState.
                eng_offdes = eng

                #----- fixed parameters (tasopt-j9l.53: read design flow areas from typed DesignState)
                A2  = eng_offdes.design.A2
                A2_5 = eng_offdes.design.A25
                A5  = eng_offdes.design.A8
                A7  = eng_offdes.design.A18

                # tasopt-j9l.23: read map anchors from typed DesignState
                Nb_fan_des  = eng_offdes.design.Nb_fan_des
                Nb_lpc_des = eng_offdes.design.Nb_lpc_des
                Nb_hpc_des = eng_offdes.design.Nb_hpc_des
                Nb_hpt_des = eng_offdes.design.Nb_hpt_des
                Nb_lpt_des = eng_offdes.design.Nb_lpt_des

                mb_fan_des  = eng_offdes.design.mb_fan_des
                mb_lpc_des = eng_offdes.design.mb_lpc_des
                mb_hpc_des = eng_offdes.design.mb_hpc_des
                mb_hpt_des = eng_offdes.design.mb_hpt_des
                mb_lpt_des = eng_offdes.design.mb_lpt_des

                pi_fan_des  = eng_offdes.design.pi_fan_des
                pi_lpc_des = eng_offdes.design.pi_lpc_des
                pi_hpc_des = eng_offdes.design.pi_hpc_des
                pi_hpt_des = eng_offdes.design.pi_hpt_des
                pi_lpt_des = eng_offdes.design.pi_lpt_des

                if (initializes_engine)
                        #------ force TFOPER to initialize these state variables
                        mbf = 0.0
                        mblc = 0.0
                        mbhc = 0.0
                        pif = 0.0
                        pilc = 0.0
                        pihc = 0.0
                        pt8 = 0.0
                        M2 = 1.0
                        M2_5 = 1.0
                else
                        #------ use existing state variables as initial guesses
                        mbf  = eng.mbf
                        mblc = eng.mblc
                        mbhc = eng.mbhc
                        pif  = max(eng.pif,  1.1)
                        pilc = max(eng.pilc, 1.1)
                        pihc = max(eng.pihc, 1.1)
                        pt8  = eng.st8.pt
                        M2   = eng.design.M2
                        M2_5  = eng.design.M25
                end

                Fe = 0.0
                Tt4 = eng.Tt4

                if opt_calc_call == CalcMode.FixedTt4OffDes
                        #------ specified Tt4 -- Fe will be computed
                        nothing; #nothing special is done
                elseif opt_calc_call == CalcMode.FixedFeOffDes
                        #------ specified Fe -- Tt4 will be computed (set initial guess here)
                        Fe = eng.Fe
                end

                if (Lprint)
                        println(cplab[ip], opt_calc_call, Tt4, Fe)
                        println("Calling TFOPER...", opt_calc_call, opt_cooling, ip)
                        println(DAwsurf, para[iagamV])
                        println(rho0, u0, parg[igWMTO])
                        println(mcore, M2, M2_5)
                        println("Phiinl, Kinl", Phiinl, Kinl)

                end

                TSFC, Fsp, hfuel, ff,
                Fe, mcore,
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
                T2_5, u2_5, p2_5, cp2_5, R2_5, M2_5,
                T5, u5, p5, cp5, R5, M5,
                T6, u6, p6, cp6, R6, M6, A6,
                T7, u7, p7, cp7, R7, M7,
                T8, u8, p8, cp8, R8, M8, A8,
                u9, A9,
                epf, eplc, ephc, epht, eplt,
                etaf, etalc, etahc, etaht, etalt,
                Lconv = tfoper!(gee, M0, T0, p0, a0, Tref, pref,
                        Phiinl, Kinl, eng_has_BLI_cores,
                        pid, pib, pifn, pitn,
                        Gearf,
                        pi_fan_des, pi_lpc_des, pi_hpc_des, pi_hpt_des, pi_lpt_des,
                        mb_fan_des, mb_lpc_des, mb_hpc_des, mb_hpt_des, mb_lpt_des,
                        Nb_fan_des, Nb_lpc_des, Nb_hpc_des, Nb_hpt_des, Nb_lpt_des,
                        A2, A2_5, A5, A7,
                        opt_calc_call,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt25off, pt25off,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Fe,
                        M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt8, mcore, M2_5, 
                        Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                        Δp_PreC, Δp_InterC, Δp_Regen)

                # ── EXIT: populate EngineState with off-design results ────────
                # Mirror the sizing exit block (lines 290–385).  Station 19c
                # (PreCoolerOut) and 25c (InterCoolerOut) are NOT in pare but
                # ARE returned by tfoper! — the typed state captures them here.
                # For st4: Tt4 is INPUT for FixedTt4OffDes and OUTPUT for
                # FixedFeOffDes; capture the current value in both cases.
                eng_offdes.st0.Tt  = Tt0;   eng_offdes.st0.ht  = ht0
                eng_offdes.st0.pt  = pt0;   eng_offdes.st0.cpt = cpt0;   eng_offdes.st0.Rt = Rt0
                eng_offdes.st0.u   = u0

                eng_offdes.st12.Tt  = Tt12; eng_offdes.st12.ht  = ht12
                eng_offdes.st12.pt  = pt12; eng_offdes.st12.cpt = cpt12; eng_offdes.st12.Rt = Rt12

                eng_offdes.st2a.Tt  = Tt2a; eng_offdes.st2a.ht  = ht2a
                eng_offdes.st2a.pt  = pt2a; eng_offdes.st2a.cpt = cpt2a; eng_offdes.st2a.Rt = Rt2a

                # st2ac (PreCoolerOut) — returned by tfoper! but NOT in pare
                eng_offdes.st2ac.Tt  = Tt2ac; eng_offdes.st2ac.ht  = ht2ac
                eng_offdes.st2ac.pt  = pt2ac; eng_offdes.st2ac.cpt = cpt2ac; eng_offdes.st2ac.Rt = Rt2ac

                eng_offdes.st2.Tt  = Tt2;  eng_offdes.st2.ht  = ht2
                eng_offdes.st2.pt  = pt2;  eng_offdes.st2.cpt = cpt2;  eng_offdes.st2.Rt = Rt2
                eng_offdes.st2.Ts  = T2;   eng_offdes.st2.ps  = p2
                eng_offdes.st2.Rs  = R2;   eng_offdes.st2.cps = cp2;   eng_offdes.st2.u  = u2
                eng_offdes.st2.A   = A2;   eng_offdes.st2.mdot = mcore

                eng_offdes.st13.Tt  = Tt13; eng_offdes.st13.ht  = ht13
                eng_offdes.st13.pt  = pt13; eng_offdes.st13.cpt = cpt13; eng_offdes.st13.Rt = Rt13

                eng_offdes.st25.Tt  = Tt2_5; eng_offdes.st25.ht  = ht2_5
                eng_offdes.st25.pt  = pt2_5; eng_offdes.st25.cpt = cpt2_5; eng_offdes.st25.Rt = Rt2_5
                eng_offdes.st25.Ts  = T2_5;  eng_offdes.st25.ps  = p2_5
                eng_offdes.st25.Rs  = R2_5;  eng_offdes.st25.cps = cp2_5;  eng_offdes.st25.u  = u2_5
                eng_offdes.st25.A   = A2_5

                # st2_5c (InterCoolerOut) — returned by tfoper! but NOT in pare
                eng_offdes.st25c.Tt  = Tt2_5c; eng_offdes.st25c.ht  = ht2_5c
                eng_offdes.st25c.pt  = pt2_5c; eng_offdes.st25c.cpt = cpt2_5c; eng_offdes.st25c.Rt = Rt2_5c

                eng_offdes.st3.Tt  = Tt3;  eng_offdes.st3.ht  = ht3
                eng_offdes.st3.pt  = pt3;  eng_offdes.st3.cpt = cpt3;  eng_offdes.st3.Rt = Rt3

                eng_offdes.st4.Tt  = Tt4;  eng_offdes.st4.ht  = ht4
                eng_offdes.st4.pt  = pt4;  eng_offdes.st4.cpt = cpt4;  eng_offdes.st4.Rt = Rt4

                eng_offdes.st41.Tt  = Tt4_1; eng_offdes.st41.ht  = ht4_1
                eng_offdes.st41.pt  = pt4_1; eng_offdes.st41.cpt = cpt4_1; eng_offdes.st41.Rt = Rt4_1

                eng_offdes.st45.Tt  = Tt4_5; eng_offdes.st45.ht  = ht4_5
                eng_offdes.st45.pt  = pt4_5; eng_offdes.st45.cpt = cpt4_5; eng_offdes.st45.Rt = Rt4_5

                eng_offdes.st5.Tt  = Tt5; eng_offdes.st5.ht  = ht5
                eng_offdes.st5.pt  = pt5; eng_offdes.st5.cpt = cpt5; eng_offdes.st5.Rt = Rt5

                eng_offdes.st8.Tt  = Tt8;  eng_offdes.st8.ht  = ht8
                eng_offdes.st8.pt  = pt8;  eng_offdes.st8.cpt = cpt8;  eng_offdes.st8.Rt = Rt8
                eng_offdes.st8.Ts  = T5;   eng_offdes.st8.ps  = p5
                eng_offdes.st8.Rs  = R5;   eng_offdes.st8.cps = cp5;   eng_offdes.st8.u  = u5
                eng_offdes.st8.A   = A5

                eng_offdes.st9.Ts  = T6;   eng_offdes.st9.ps  = p6
                eng_offdes.st9.Rs  = R6;   eng_offdes.st9.cps = cp6;   eng_offdes.st9.u  = u6
                eng_offdes.st9.A   = A6

                eng_offdes.st18.Tt  = Tt18;  eng_offdes.st18.ht  = ht18
                eng_offdes.st18.pt  = pt18;  eng_offdes.st18.cpt = cpt18;  eng_offdes.st18.Rt = Rt18
                eng_offdes.st18.Ts  = T7;   eng_offdes.st18.ps  = p7
                eng_offdes.st18.Rs  = R7;   eng_offdes.st18.cps = cp7;   eng_offdes.st18.u  = u7
                eng_offdes.st18.A   = A7

                eng_offdes.st19.Ts  = T8;   eng_offdes.st19.ps  = p8
                eng_offdes.st19.Rs  = R8;   eng_offdes.st19.cps = cp8;   eng_offdes.st19.u  = u8
                eng_offdes.st19.A   = A8

                eng_offdes.st25off.u   = u9;   eng_offdes.st25off.A   = A9

                # Cooling — mirror the sizing EXIT block (lines 376-378) for both
                # operating modes (FixedCoolingFlowRatio and FixedTmetal).
                # epsrow/Tmrow are populated by tfoper! in place; mofft/mcore
                # are available as tfoper! return values.
                eng_offdes.design.epsrow = SVector{4,Float64}(epsrow[1], epsrow[2], epsrow[3], epsrow[4])
                eng_offdes.design.Tmrow  = SVector{4,Float64}(Tmrow[1],  Tmrow[2],  Tmrow[3],  Tmrow[4])
                eng_offdes.design.fc     = (1.0 - mofft / mcore) * sum(epsrow)

                # Performance rollup scalars (tasopt-j9l.52)
                # BPR is computed from off-design corrected flows (not the design-point input).
                # Fe is the per-engine net thrust (OUTPUT in FixedTt4; INPUT retained in FixedFe).
                eng_offdes.BPR   = mbf / mblc * sqrt(Tt2ac / Tt2) * pt2 / pt2ac
                eng_offdes.Fe    = Fe
                eng_offdes.TSFC  = TSFC
                eng_offdes.Fsp   = Fsp
                eng_offdes.mfuel = ff * mcore * neng

                # Component adiabatic efficiencies (tasopt-j9l.63.1)
                eng_offdes.etaf  = etaf
                eng_offdes.etalc = etalc
                eng_offdes.etahc = etahc
                eng_offdes.etaht = etaht
                eng_offdes.etalt = etalt

                # Overall propulsion efficiencies (tasopt-j9l.63.2)
                # Guard: zero at ground-idle (Fe ≤ 0 or Fsp ≤ 0).
                if Fe > 0.0 && Fsp > 0.0
                    _eta_overall = Fe * u0 / (ff * mcore * hfuel)
                    _eta_prop    = 2.0 / (2.0 + Fsp)
                    eng_offdes.eta_overall = _eta_overall
                    eng_offdes.eta_prop    = _eta_prop
                    eng_offdes.eta_thermal = _eta_overall / _eta_prop
                else
                    eng_offdes.eta_overall = 0.0
                    eng_offdes.eta_prop    = 0.0
                    eng_offdes.eta_thermal = 0.0
                end

                # Scalar fields written directly into typed EngineState (tasopt-j9l.45.14.1).
                eng_offdes.hfuel  = hfuel
                eng_offdes.ff     = ff
                eng_offdes.mofft  = mofft
                eng_offdes.Pofft  = Pofft
                eng_offdes.Phiinl = Phiinl
                eng_offdes.Kinl   = Kinl
                eng_offdes.Nf    = Nbf  * sqrt(Tt2   / Tref)
                eng_offdes.N1    = Nblc * sqrt(Tt2ac / Tref)
                eng_offdes.N2    = Nbhc * sqrt(Tt2_5c / Tref)
                eng_offdes.Nbf   = Nbf;   eng_offdes.Nblc = Nblc;  eng_offdes.Nbhc = Nbhc
                eng_offdes.mbf   = mbf;   eng_offdes.mblc = mblc;  eng_offdes.mbhc = mbhc
                eng_offdes.pif   = pif;   eng_offdes.pilc = pilc;  eng_offdes.pihc = pihc
                eng_offdes.epf   = epf;   eng_offdes.eplc = eplc;  eng_offdes.ephc = ephc
                eng_offdes.epht  = epht;  eng_offdes.eplt = eplt
                # M2/M2_5: off-design converged values used as initial guess for next call.
                eng_offdes.design.M2  = M2
                eng_offdes.design.M25 = M2_5

                if (Lprint)
                        println("exited TFOPER", Lconv)
                end

                if (!Lconv)
                        #@warn "Convergence failed on operating point: $ip"
                        eng_offdes.ConvFail = 1.0
                end

                fo = mofft / mcore

                # println("exited TFOPER call")

        end
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        # Cooling (tasopt-j9l.20):
        # FixedCoolingFlowRatio — Tmrow is written into eng.design.Tmrow in EXIT block.
        # FixedCoolingFlowRatio — epsrow is an INPUT (read from eng.design.epsrow at entry).
        # FixedTmetal — epsrow is written into eng.design.epsrow in EXIT block.
        # FixedTmetal — fc (iefc) is a COMPUTED OUTPUT only in this mode; written via
        #   eng.design.fc (set in EXIT block at line 656).
        # bare-pare outputs removed (tasopt-j9l.45.14.1/52/63.1): written directly into typed EngineState.

        if (M5 <= 0.999999)
                ichoke5 = 0
        else
                ichoke5 = 1
        end

        if (M7 <= 0.999999)
                ichoke7 = 0
        else
                ichoke7 = 1
        end

        return ichoke5, ichoke7
end # tfcalc

"""
    _update_engine_sizing!(eng, result, BPR, Fe, mofft, Pofft, Phiinl, Kinl, neng,
                           pif, pilc, pihc, ruc, M4a, epolht, epollt, Gearf, parg)

Write the outputs of `tfsize!` (collected in `result::SizingResult`) back into the
typed `EngineState` `eng`.  Also computes the derived design-point map anchors
(corrected mass flows, corrected spool speeds, turbine pressure ratios) and writes
the fan/compressor geometric diameters to `parg`.

This helper centralises all station and scalar writeback for the `CalcMode.Sizing`
branch of `tfcalc!`.  Every field mapping is documented inline so that adding a new
station or scalar requires editing only this function.

Mapping contract (station fields → EngineState):
  result.Tt0/ht0/pt0/cpt0/Rt0/u0  → eng.st0
  result.Tt12..Rt12               → eng.st12
  result.Tt2a..Rt2a               → eng.st2a
  result.Tt2ac..Rt2ac             → eng.st2ac  (PreCoolerOut; not in bare-pare)
  result.Tt2..Rt2 + T2..A2        → eng.st2 (totals + statics + area + mdot)
  result.Tt13..Rt13               → eng.st13
  result.Tt2_5..Rt2_5 + T2_5..A2_5 → eng.st25 (totals + statics + area)
  result.Tt2_5c..Rt2_5c           → eng.st25c (InterCoolerOut; not in bare-pare)
  result.Tt3..Rt3                 → eng.st3
  result.ht4/pt4/cpt4/Rt4         → eng.st4 (Tt4 is an input; not overwritten)
  result.Tt4_1..Rt4_1             → eng.st41
  result.Tt4_5..Rt4_5             → eng.st45
  result.Tt5..Rt5                 → eng.st5
  result.Tt8..Rt8 + T5..A5        → eng.st8 (totals + statics mapped from LPT exit)
  result.T6..A6                   → eng.st9 (statics only)
  result.Tt18..Rt18 + T7..A7      → eng.st18 (totals + statics)
  result.T8..A8                   → eng.st19 (statics only)
  result.u9/A9                    → eng.st25off
"""
function _update_engine_sizing!(eng::EngineState, result::SizingResult,
        BPR::Float64, Fe::Float64, mofft::Float64, Pofft::Float64,
        Phiinl::Float64, Kinl::Float64, neng::Float64,
        pif::Float64, pilc::Float64, pihc::Float64,
        ruc::Float64, M4a::Float64, epolht::Float64, epollt::Float64,
        Gearf::Float64, parg::Vector{Float64})

        fo = mofft / result.mcore

        # Corrected mass flows
        mbf  = result.mcore * sqrt(result.Tt2   / Tref) / (result.pt2   / pref) * BPR
        mblc = result.mcore * sqrt(result.Tt2ac / Tref) / (result.pt2ac / pref)
        mbhc = result.mcore * sqrt(result.Tt2_5c / Tref) / (result.pt2_5c / pref) * (1.0 - fo)
        mbht = result.mcore * sqrt(result.Tt4_1 / Tref) / (result.pt4_1 / pref) * (1.0 - fo + result.ff)
        mblt = result.mcore * sqrt(result.Tt4_5 / Tref) / (result.pt4_5 / pref) * (1.0 - fo + result.ff)

        # Spool speed fractions: design point = 1 by definition
        Nf = 1.0 / Gearf
        N1 = 1.0
        N2 = 1.0

        Nbf  = Nf / sqrt(result.Tt2   / Tref)
        Nblc = N1 / sqrt(result.Tt2ac / Tref)
        Nbhc = N2 / sqrt(result.Tt2_5c / Tref)
        Nbht = N2 / sqrt(result.Tt4_1 / Tref)
        Nblt = N1 / sqrt(result.Tt4_5 / Tref)

        # Map anchors (design-point corrected quantities)
        mb_fan_des = mbf;  mb_lpc_des = mblc;  mb_hpc_des = mbhc
        mb_hpt_des = mbht; mb_lpt_des = mblt

        pi_fan_des = pif;  pi_lpc_des = pilc;  pi_hpc_des = pihc
        pi_hpt_des = result.pt4_1 / result.pt4_5
        pi_lpt_des = result.pt4_5 / result.pt5

        Nb_fan_des = Nbf;  Nb_lpc_des = Nblc;  Nb_hpc_des = Nbhc
        Nb_hpt_des = Nbht; Nb_lpt_des = Nblt

        # Recalculate turbine pressure ratios using the slightly approximate form
        # consistent with tfoper!'s turbine efficiency function.
        Trh = result.Tt4_1 / (result.Tt4_1 + (result.ht4_5 - result.ht4_1) / result.cpt4_1)
        Trl = result.Tt4_5 / (result.Tt4_5 + (result.ht5   - result.ht4_5) / result.cpt4_5)
        gexh = result.cpt4_1 / (result.Rt4_1 * epolht)
        gexl = result.cpt4_5 / (result.Rt4_5 * epollt)
        pi_hpt_des = Trh^gexh
        pi_lpt_des = Trl^gexl

        # ── Station total and static states ──────────────────────────────────
        eng.st0.Tt  = result.Tt0;  eng.st0.ht  = result.ht0
        eng.st0.pt  = result.pt0;  eng.st0.cpt = result.cpt0;  eng.st0.Rt = result.Rt0
        eng.st0.u   = result.u0

        eng.st12.Tt  = result.Tt12; eng.st12.ht  = result.ht12
        eng.st12.pt  = result.pt12; eng.st12.cpt = result.cpt12; eng.st12.Rt = result.Rt12

        eng.st2a.Tt  = result.Tt2a; eng.st2a.ht  = result.ht2a
        eng.st2a.pt  = result.pt2a; eng.st2a.cpt = result.cpt2a; eng.st2a.Rt = result.Rt2a

        # st2ac (PreCoolerOut) — returned by tfsize! but not in bare-pare
        eng.st2ac.Tt  = result.Tt2ac; eng.st2ac.ht  = result.ht2ac
        eng.st2ac.pt  = result.pt2ac; eng.st2ac.cpt = result.cpt2ac; eng.st2ac.Rt = result.Rt2ac

        eng.st2.Tt  = result.Tt2;  eng.st2.ht  = result.ht2
        eng.st2.pt  = result.pt2;  eng.st2.cpt = result.cpt2;  eng.st2.Rt = result.Rt2
        eng.st2.Ts  = result.T2;   eng.st2.ps  = result.p2
        eng.st2.Rs  = result.R2;   eng.st2.cps = result.cp2;   eng.st2.u = result.u2
        eng.st2.A   = result.A2;   eng.st2.mdot = result.mcore

        eng.st13.Tt  = result.Tt13; eng.st13.ht  = result.ht13
        eng.st13.pt  = result.pt13; eng.st13.cpt = result.cpt13; eng.st13.Rt = result.Rt13

        eng.st25.Tt  = result.Tt2_5; eng.st25.ht  = result.ht2_5
        eng.st25.pt  = result.pt2_5; eng.st25.cpt = result.cpt2_5; eng.st25.Rt = result.Rt2_5
        eng.st25.Ts  = result.T2_5;  eng.st25.ps  = result.p2_5
        eng.st25.Rs  = result.R2_5;  eng.st25.cps = result.cp2_5;  eng.st25.u = result.u2_5
        eng.st25.A   = result.A2_5

        # st25c (InterCoolerOut) — returned by tfsize! but not in bare-pare
        eng.st25c.Tt  = result.Tt2_5c; eng.st25c.ht  = result.ht2_5c
        eng.st25c.pt  = result.pt2_5c; eng.st25c.cpt = result.cpt2_5c; eng.st25c.Rt = result.Rt2_5c

        eng.st3.Tt  = result.Tt3;  eng.st3.ht  = result.ht3
        eng.st3.pt  = result.pt3;  eng.st3.cpt = result.cpt3;  eng.st3.Rt = result.Rt3

        # st4: Tt4 is an INPUT (not overwritten); ht4/pt4/cpt4/Rt4 are sizing outputs
        eng.st4.ht  = result.ht4;  eng.st4.pt  = result.pt4
        eng.st4.cpt = result.cpt4; eng.st4.Rt  = result.Rt4

        eng.st41.Tt  = result.Tt4_1; eng.st41.ht  = result.ht4_1
        eng.st41.pt  = result.pt4_1; eng.st41.cpt = result.cpt4_1; eng.st41.Rt = result.Rt4_1

        eng.st45.Tt  = result.Tt4_5; eng.st45.ht  = result.ht4_5
        eng.st45.pt  = result.pt4_5; eng.st45.cpt = result.cpt4_5; eng.st45.Rt = result.Rt4_5

        eng.st5.Tt  = result.Tt5; eng.st5.ht  = result.ht5
        eng.st5.pt  = result.pt5; eng.st5.cpt = result.cpt5; eng.st5.Rt = result.Rt5

        # st8 totals from station 8; statics from LPT exit (station 5) — existing convention
        eng.st8.Tt  = result.Tt8;  eng.st8.ht  = result.ht8
        eng.st8.pt  = result.pt8;  eng.st8.cpt = result.cpt8;  eng.st8.Rt = result.Rt8
        eng.st8.Ts  = result.T5;   eng.st8.ps  = result.p5
        eng.st8.Rs  = result.R5;   eng.st8.cps = result.cp5;   eng.st8.u = result.u5
        eng.st8.A   = result.A5

        eng.st9.Ts  = result.T6;   eng.st9.ps  = result.p6
        eng.st9.Rs  = result.R6;   eng.st9.cps = result.cp6;   eng.st9.u = result.u6
        eng.st9.A   = result.A6

        eng.st18.Tt  = result.Tt18;  eng.st18.ht  = result.ht18
        eng.st18.pt  = result.pt18;  eng.st18.cpt = result.cpt18;  eng.st18.Rt = result.Rt18
        eng.st18.Ts  = result.T7;   eng.st18.ps  = result.p7
        eng.st18.Rs  = result.R7;   eng.st18.cps = result.cp7;   eng.st18.u = result.u7
        eng.st18.A   = result.A7

        eng.st19.Ts  = result.T8;   eng.st19.ps  = result.p8
        eng.st19.Rs  = result.R8;   eng.st19.cps = result.cp8;   eng.st19.u = result.u8
        eng.st19.A   = result.A8

        eng.st25off.u = result.u9; eng.st25off.A = result.A9

        # ── DesignState — frozen scalars for all off-design calls ────────────
        eng.design.pi_fan_des  = pi_fan_des;  eng.design.pi_lpc_des = pi_lpc_des
        eng.design.pi_hpc_des  = pi_hpc_des;  eng.design.pi_hpt_des = pi_hpt_des
        eng.design.pi_lpt_des  = pi_lpt_des
        eng.design.mb_fan_des  = mb_fan_des;  eng.design.mb_lpc_des = mb_lpc_des
        eng.design.mb_hpc_des  = mb_hpc_des;  eng.design.mb_hpt_des = mb_hpt_des
        eng.design.mb_lpt_des  = mb_lpt_des
        eng.design.Nb_fan_des  = Nb_fan_des;  eng.design.Nb_lpc_des = Nb_lpc_des
        eng.design.Nb_hpc_des  = Nb_hpc_des;  eng.design.Nb_hpt_des = Nb_hpt_des
        eng.design.Nb_lpt_des  = Nb_lpt_des
        eng.design.A2   = result.A2;  eng.design.A25  = result.A2_5
        eng.design.A8   = result.A5;  eng.design.A18  = result.A7
        eng.design.epsrow = SVector{4,Float64}(result.epsrow[1], result.epsrow[2],
                                               result.epsrow[3], result.epsrow[4])
        eng.design.Tmrow  = SVector{4,Float64}(result.Tmrow[1],  result.Tmrow[2],
                                               result.Tmrow[3],  result.Tmrow[4])
        eng.design.fc   = (1.0 - fo) * sum(result.epsrow)
        eng.design.ruc  = ruc
        eng.design.M4a  = M4a

        # ── Performance rollup scalars ────────────────────────────────────────
        # BPR and Fe are INPUTs to tfsize! (passed through unchanged).
        eng.TSFC  = result.TSFC
        eng.Fe    = Fe
        eng.Fsp   = result.Fsp
        eng.BPR   = BPR
        eng.mfuel = result.ff * result.mcore * neng

        # Component adiabatic efficiencies (tasopt-j9l.63.1)
        eng.etaf  = result.etaf
        eng.etalc = result.etalc
        eng.etahc = result.etahc
        eng.etaht = result.etaht
        eng.etalt = result.etalt

        # Overall propulsion efficiencies (tasopt-j9l.63.2)
        # Guard: zero at ground-idle (Fe ≤ 0 or Fsp ≤ 0).
        if Fe > 0.0 && result.Fsp > 0.0
                _eta_overall = Fe * result.u0 / (result.ff * result.mcore * result.hfuel)
                _eta_prop    = 2.0 / (2.0 + result.Fsp)
                eng.eta_overall = _eta_overall
                eng.eta_prop    = _eta_prop
                eng.eta_thermal = _eta_overall / _eta_prop
        else
                eng.eta_overall = 0.0
                eng.eta_prop    = 0.0
                eng.eta_thermal = 0.0
        end

        # Scalar operating outputs
        eng.hfuel  = result.hfuel
        eng.ff     = result.ff
        eng.mofft  = mofft
        eng.Pofft  = Pofft
        eng.Phiinl = Phiinl
        eng.Kinl   = Kinl
        eng.Nf    = Nf;    eng.N1   = N1;    eng.N2   = N2
        eng.Nbf   = Nbf;   eng.Nblc = Nblc;  eng.Nbhc = Nbhc
        eng.mbf   = mbf;   eng.mblc = mblc;  eng.mbhc = mbhc
        eng.pif   = pif;   eng.pilc = pilc;  eng.pihc = pihc
        eng.epf   = result.epf;  eng.eplc = result.eplc;  eng.ephc = result.ephc
        eng.epht  = result.epht; eng.eplt = result.eplt

        # Fan/compressor diameters written back to parg
        HTRf  = parg[igHTRf]
        HTRlc = parg[igHTRlc]
        HTRhc = parg[igHTRhc]
        Alc = result.A2 / (1.0 + BPR)
        parg[igdfan]   = sqrt(4.0 * result.A2  / (pi * (1.0 - HTRf^2)))
        parg[igdlcomp] = sqrt(4.0 * Alc         / (pi * (1.0 - HTRlc^2)))
        parg[igdhcomp] = sqrt(4.0 * result.A2_5 / (pi * (1.0 - HTRhc^2)))
        nothing
end

"""
    check_engine_convergence_failure(ac, imission)

Check whether any engine operating point failed to converge for `imission`.
Reads `ConvFail` from the typed `EngineState` of each mission point.
Emits a warning if any point did not converge.
"""
function check_engine_convergence_failure(ac, imission::Int)
        points = ac.missions[imission].points
        failed = [ip for (ip, pt) in enumerate(points) if pt.engine.ConvFail != 0.0]
        if !isempty(failed)
                @warn "Some engine points failed to converge: point indices $failed"
        end
end
