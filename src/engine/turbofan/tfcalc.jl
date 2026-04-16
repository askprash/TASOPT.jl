""" 
    tfcalc!(wing, engine, parg, para, pare, ip, ifuel, opt_calc_call, opt_cooling, initializes_engine)

Calls on-design sizing function [`tfsize!`](@ref) or off-design analysis function
[`tfoper!`](@ref) for one operating point `ip`.

!!! details "🔃 Inputs and Outputs"
    **Input:**
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
function tfcalc!(wing, engine, parg::Vector{Float64}, para, pare, ip::Int64, ifuel::Int64,
        opt_calc_call::CalcMode.T, opt_cooling::CoolingOpt.T, initializes_engine::Bool)

        # ── ENTRY: build typed EngineState from pare once, shared by both paths ─
        # Ambient scalars and all station thermodynamics are populated here so
        # that the sizing and off-design branches can alias `eng` rather than
        # each constructing their own copy.
        eng = EngineState{Float64}()
        pare_to_engine_state!(eng, pare)

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

        TSFC = pare[ieTSFC]
        Fsp = pare[ieFsp]
        hfuel = pare[iehfuel]
        Tfuel = pare[ieTfuel]
        Tt4 = pare[ieTt4]
        BPR = pare[ieBPR]
        pif = pare[iepif]
        pilc = pare[iepilc]
        pihc = pare[iepihc]
        pid = pare[iepid]
        pib = pare[iepib]
        pifn = pare[iepifn]
        pitn = pare[iepitn]
        epolf = pare[ieepolf]
        epollc = pare[ieepollc]
        epolhc = pare[ieepolhc]
        epolht = pare[ieepolht]
        epollt = pare[ieepollt]
        etab = pare[ieetab]
        M2 = pare[ieM2]
        M25 = pare[ieM25]
        # Ambient / freestream scalars — read from typed state (not pare).
        # Tt0/ht0/pt0/cpt0/Rt0 are computed outputs of tfsize!/tfoper! and are
        # not pre-read here; they are assigned from the returned tuple below.
        M0   = eng.M0
        T0   = eng.T0
        p0   = eng.p0
        a0   = eng.a0
        u0   = eng.st0.u
        rho0 = pare[ierho0]   # not in EngineState; used only for BLI
        mu0  = pare[iemu0]    # not in EngineState
        mcore = eng.st2.mdot
        dTstrk = pare[iedTstrk]
        StA = pare[ieStA]
        Mtexit = pare[ieMtexit]
        M4a = pare[ieM4a]
        ruc = pare[ieruc]
        efilm = pare[ieefilm]
        tfilm = pare[ietfilm]
        epsl = pare[ieepsl]
        epsh = pare[ieepsh]
        epsrow = zeros(ncrowx)
        Tmrow = zeros(ncrowx)
        hvap = pare[iehvapcombustor] #Enthalpy of vaporization of fuel

        #Effect of cooling on HPT efficiency
        epht_fc = pare[iedehtdfc]
        fc0 = fc0 = pare[iefc0]

        #Heat exchanger variables
        Δh_PreC = pare[iePreCDeltah]
        Δh_InterC = pare[ieInterCDeltah]
        Δh_Regen = pare[ieRegenDeltah]
        Δh_TurbC = pare[ieTurbCDeltah]
        Δp_PreC = pare[iePreCDeltap]
        Δp_InterC = pare[ieInterCDeltap]
        Δp_Regen = pare[ieRegenDeltap]

        if opt_cooling == CoolingOpt.FixedCoolingFlowRatio
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        epsrow[icrow] = eng.design.epsrow[icrow]
                end
        elseif opt_cooling == CoolingOpt.FixedTmetal
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        #cc      Tmrow[icrow]  = pare(ieTmet1+icrow-1)
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
        Pofft_HX = pare[ieHXrecircP] #power offtakes to drive heat exchanger recirculation per engine
        mofft = (mofWpay * Wpay + mofWMTO * WMTO) / neng
        Pofft = (PofWpay * Wpay + PofWMTO * WMTO) / neng + Pofft_HX

        Tt9 = pare[ieTt9]
        pt9 = pare[iept9]

        #--------------------------------------------------------------------------
        #Engine model convergence
        pare[ieConvFail] = 0.0 #Converged by default

        # #--------------------------------------------------------------------------
        if opt_calc_call == CalcMode.Sizing
                #----- engine sizing case

                # ── ENTRY ────────────────────────────────────────────────────
                # Alias the top-level EngineState built before the if/else.
                # After tfsize! returns we populate every station and the
                # DesignState from local variables, then project back to pare.
                eng_design = eng

                Fe = pare[ieFe]

                if (Lprint)
                        println("TFSIZE  M0 p0 =", M0, p0)
                        println("     pif pilc =", pif, pilc)
                        println("        Tt4 F =", Tt4, Fe)
                        println("ncrow =", ncrow)
                        println("        alt_m =", para[iaalt])
                end

                epsrow, Tmrow,
                TSFC, Fsp, hfuel, ff, mcore,
                Tt0, ht0, pt0, cpt0, Rt0,
                Tt18, ht18, pt18, cpt18, Rt18,
                Tt19, ht19, pt19, cpt19, Rt19,
                Tt19c, ht19c, pt19c, cpt19c, Rt19c,
                Tt2, ht2, pt2, cpt2, Rt2,
                Tt21, ht21, pt21, cpt21, Rt21,
                Tt25, ht25, pt25, cpt25, Rt25,
                Tt25c, ht25c, pt25c, cpt25c, Rt25c,
                Tt3, ht3, pt3, cpt3, Rt3,
                ht4, pt4, cpt4, Rt4,
                Tt41, ht41, pt41, cpt41, Rt41,
                Tt45, ht45, pt45, cpt45, Rt45,
                Tt49, ht49, pt49, cpt49, Rt49,
                Tt5, ht5, pt5, cpt5, Rt5,
                Tt7, ht7, pt7, cpt7, Rt7,
                u0,
                T2, u2, p2, cp2, R2, A2,
                T25, u25, p25, cp25, R25, A25,
                T5, u5, p5, cp5, R5, A5,
                T6, u6, p6, cp6, R6, A6,
                T7, u7, p7, cp7, R7, A7,
                T8, u8, p8, cp8, R8, A8,
                u9, A9,
                epf, eplc, ephc, epht, eplt,
                etaf, etalc, etahc, etaht, etalt,
                Lconv = tfsize!(gee, M0, T0, p0, a0, M2, M25,
                        Fe, Phiinl, Kinl, eng_has_BLI_cores,
                        BPR, pif, pilc, pihc,
                        pid, pib, pifn, pitn,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt9, pt9, Tt4,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                        Δp_PreC, Δp_InterC, Δp_Regen)

                #        tadd(time0,t_tfsize)

                M5 = u5 / sqrt(T5 * R5 * cp5 / (cp5 - R5))
                M6 = u6 / sqrt(T6 * R6 * cp6 / (cp6 - R6))
                M7 = u7 / sqrt(T7 * R7 * cp7 / (cp7 - R7))
                M8 = u8 / sqrt(T8 * R8 * cp8 / (cp8 - R8))

                #      pare(ieM5 ) = M5   
                #      pare(ieM6 ) = M6   
                #      pare(ieM7 ) = M7   
                #      pare(ieM8 ) = M8   

                fo = mofft / mcore

                if (Lprint)
                        println("exited TFSIZE", fo)
                end

                #----- corrected mass flows
                mbf = mcore * sqrt(Tt2 / Tref) / (pt2 / pref) * BPR
                mblc = mcore * sqrt(Tt19c / Tref) / (pt19c / pref)
                mbhc = mcore * sqrt(Tt25c / Tref) / (pt25c / pref) * (1.0 - fo)
                mbht = mcore * sqrt(Tt41 / Tref) / (pt41 / pref) * (1.0 - fo + ff)
                mblt = mcore * sqrt(Tt45 / Tref) / (pt45 / pref) * (1.0 - fo + ff)

                #----- spool speed fractions for design case are unity by definition
                Nf = 1.0 / Gearf
                N1 = 1.0
                N2 = 1.0

                Nbf = Nf / sqrt(Tt2 / Tref)
                Nblc = N1 / sqrt(Tt19c / Tref)
                Nbhc = N2 / sqrt(Tt25c / Tref)
                Nbht = N2 / sqrt(Tt41 / Tref)
                Nblt = N1 / sqrt(Tt45 / Tref)

                #----- set quantities fixed by this design case for all operating points
                mbfD = mbf
                mblcD = mblc
                mbhcD = mbhc
                mbhtD = mbht
                mbltD = mblt

                pifD = pif
                pilcD = pilc
                pihcD = pihc
                pihtD = pt41 / pt45
                piltD = pt45 / pt49

                NbfD = Nbf
                NblcD = Nblc
                NbhcD = Nbhc
                NbhtD = Nbht
                NbltD = Nblt

                #----- but recalculate turbine pressure ratios using slightly approximate form,
                #-      to be fuly consistent with TFOPER's turbine efficiency function
                Trh = Tt41 / (Tt41 + (ht45 - ht41) / cpt41)
                Trl = Tt45 / (Tt45 + (ht49 - ht45) / cpt45)
                gexh = cpt41 / (Rt41 * epolht)
                gexl = cpt45 / (Rt45 * epollt)
                pihtD = Trh^gexh
                piltD = Trl^gexl

                # ── EXIT: populate EngineState with sizing results ───────────
                # Station total states
                eng_design.st0.Tt  = Tt0;  eng_design.st0.ht  = ht0
                eng_design.st0.pt  = pt0;  eng_design.st0.cpt = cpt0;  eng_design.st0.Rt = Rt0
                eng_design.st0.u   = u0

                eng_design.st18.Tt  = Tt18; eng_design.st18.ht  = ht18
                eng_design.st18.pt  = pt18; eng_design.st18.cpt = cpt18; eng_design.st18.Rt = Rt18

                eng_design.st19.Tt  = Tt19; eng_design.st19.ht  = ht19
                eng_design.st19.pt  = pt19; eng_design.st19.cpt = cpt19; eng_design.st19.Rt = Rt19

                # st19c (PreCoolerOut) populated by tfsize! but not in pare — update typed state
                eng_design.st19c.Tt  = Tt19c; eng_design.st19c.ht  = ht19c
                eng_design.st19c.pt  = pt19c; eng_design.st19c.cpt = cpt19c; eng_design.st19c.Rt = Rt19c

                eng_design.st2.Tt  = Tt2;  eng_design.st2.ht  = ht2
                eng_design.st2.pt  = pt2;  eng_design.st2.cpt = cpt2;  eng_design.st2.Rt = Rt2
                eng_design.st2.Ts  = T2;   eng_design.st2.ps  = p2
                eng_design.st2.Rs  = R2;   eng_design.st2.cps = cp2;   eng_design.st2.u = u2
                eng_design.st2.A   = A2;   eng_design.st2.mdot = mcore

                eng_design.st21.Tt  = Tt21; eng_design.st21.ht  = ht21
                eng_design.st21.pt  = pt21; eng_design.st21.cpt = cpt21; eng_design.st21.Rt = Rt21

                eng_design.st25.Tt  = Tt25; eng_design.st25.ht  = ht25
                eng_design.st25.pt  = pt25; eng_design.st25.cpt = cpt25; eng_design.st25.Rt = Rt25
                eng_design.st25.Ts  = T25;  eng_design.st25.ps  = p25
                eng_design.st25.Rs  = R25;  eng_design.st25.cps = cp25;  eng_design.st25.u = u25
                eng_design.st25.A   = A25

                # st25c (InterCoolerOut) populated by tfsize! but not in pare — update typed state
                eng_design.st25c.Tt  = Tt25c; eng_design.st25c.ht  = ht25c
                eng_design.st25c.pt  = pt25c; eng_design.st25c.cpt = cpt25c; eng_design.st25c.Rt = Rt25c

                eng_design.st3.Tt  = Tt3;  eng_design.st3.ht  = ht3
                eng_design.st3.pt  = pt3;  eng_design.st3.cpt = cpt3;  eng_design.st3.Rt = Rt3

                # st4: Tt4 is an INPUT (unchanged); ht4/pt4/cpt4/Rt4 are outputs
                eng_design.st4.ht  = ht4;  eng_design.st4.pt  = pt4
                eng_design.st4.cpt = cpt4; eng_design.st4.Rt  = Rt4

                eng_design.st41.Tt  = Tt41; eng_design.st41.ht  = ht41
                eng_design.st41.pt  = pt41; eng_design.st41.cpt = cpt41; eng_design.st41.Rt = Rt41

                eng_design.st45.Tt  = Tt45; eng_design.st45.ht  = ht45
                eng_design.st45.pt  = pt45; eng_design.st45.cpt = cpt45; eng_design.st45.Rt = Rt45

                eng_design.st49.Tt  = Tt49; eng_design.st49.ht  = ht49
                eng_design.st49.pt  = pt49; eng_design.st49.cpt = cpt49; eng_design.st49.Rt = Rt49

                eng_design.st5.Tt  = Tt5;  eng_design.st5.ht  = ht5
                eng_design.st5.pt  = pt5;  eng_design.st5.cpt = cpt5;  eng_design.st5.Rt = Rt5
                eng_design.st5.Ts  = T5;   eng_design.st5.ps  = p5
                eng_design.st5.Rs  = R5;   eng_design.st5.cps = cp5;   eng_design.st5.u = u5
                eng_design.st5.A   = A5

                eng_design.st6.Ts  = T6;   eng_design.st6.ps  = p6
                eng_design.st6.Rs  = R6;   eng_design.st6.cps = cp6;   eng_design.st6.u = u6
                eng_design.st6.A   = A6

                eng_design.st7.Tt  = Tt7;  eng_design.st7.ht  = ht7
                eng_design.st7.pt  = pt7;  eng_design.st7.cpt = cpt7;  eng_design.st7.Rt = Rt7
                eng_design.st7.Ts  = T7;   eng_design.st7.ps  = p7
                eng_design.st7.Rs  = R7;   eng_design.st7.cps = cp7;   eng_design.st7.u = u7
                eng_design.st7.A   = A7

                eng_design.st8.Ts  = T8;   eng_design.st8.ps  = p8
                eng_design.st8.Rs  = R8;   eng_design.st8.cps = cp8;   eng_design.st8.u = u8
                eng_design.st8.A   = A8

                eng_design.st9.u   = u9;   eng_design.st9.A   = A9

                # DesignState — frozen scalars needed by every off-design call
                eng_design.design.pifD  = pifD;  eng_design.design.pilcD = pilcD
                eng_design.design.pihcD = pihcD; eng_design.design.pihtD = pihtD
                eng_design.design.piltD = piltD
                eng_design.design.mbfD  = mbfD;  eng_design.design.mblcD = mblcD
                eng_design.design.mbhcD = mbhcD; eng_design.design.mbhtD = mbhtD
                eng_design.design.mbltD = mbltD
                eng_design.design.NbfD  = NbfD;  eng_design.design.NblcD = NblcD
                eng_design.design.NbhcD = NbhcD; eng_design.design.NbhtD = NbhtD
                eng_design.design.NbltD = NbltD
                eng_design.design.A2    = A2;    eng_design.design.A25   = A25
                eng_design.design.A5    = A5;    eng_design.design.A7    = A7
                eng_design.design.epsrow = SVector{4,Float64}(epsrow[1], epsrow[2], epsrow[3], epsrow[4])
                eng_design.design.Tmrow  = SVector{4,Float64}(Tmrow[1],  Tmrow[2],  Tmrow[3],  Tmrow[4])
                eng_design.design.fc    = (1.0 - fo) * sum(epsrow)
                eng_design.design.ruc   = ruc
                eng_design.design.M4a   = M4a

                # Performance rollup scalars (tasopt-j9l.52)
                # BPR is an INPUT to tfsize! (unchanged from entry).
                # Fe is an INPUT to tfsize! (design thrust target).
                # TSFC/Fsp are computed outputs returned by tfsize!.
                eng_design.TSFC  = TSFC
                eng_design.Fe    = Fe
                eng_design.Fsp   = Fsp
                eng_design.BPR   = BPR
                eng_design.mfuel = ff * mcore * neng

                # ── PROJECT design scalars back to pare ──────────────────────
                # design_state_to_pare! writes frozen design-point scalars.
                # engine_state_to_pare! writes ambient scalars, station
                # thermodynamics, and performance rollup scalars; the shared
                # ambient+inlet writes below the if/else (st0/st18/st19/u0)
                # are now redundant and removed.
                design_state_to_pare!(eng_design.design, pare)
                engine_state_to_pare!(eng_design, pare)
                # pare[iemfuel] removed (tasopt-j9l.52): written by engine_state_to_pare! via eng.mfuel

                HTRf = parg[igHTRf]
                HTRlc = parg[igHTRlc]
                HTRhc = parg[igHTRhc]
                Alc = A2 / (1.0 + BPR)
                dfan = sqrt(4.0 * A2 / (pi * (1.0 - HTRf^2)))
                dlcomp = sqrt(4.0 * Alc / (pi * (1.0 - HTRlc^2)))
                dhcomp = sqrt(4.0 * A25 / (pi * (1.0 - HTRhc^2)))
                parg[igdfan] = dfan
                parg[igdlcomp] = dlcomp
                parg[igdhcomp] = dhcomp

                #--------------------------------------------------------------------------
        else
                #----- off-design operation case

                # ── ENTRY ────────────────────────────────────────────────────
                # Alias the top-level EngineState built before the if/else.
                # After tfoper! returns we populate every station from local
                # variables and project back to pare via engine_state_to_pare!.
                eng_offdes = eng

                #----- fixed parameters (tasopt-j9l.53: read design flow areas from typed DesignState)
                A2  = eng_offdes.design.A2
                A25 = eng_offdes.design.A25
                A5  = eng_offdes.design.A5
                A7  = eng_offdes.design.A7

                # tasopt-j9l.23: read map anchors from typed DesignState
                # (populated by pare_to_engine_state! at harness entry)
                NbfD  = eng_offdes.design.NbfD
                NblcD = eng_offdes.design.NblcD
                NbhcD = eng_offdes.design.NbhcD
                NbhtD = eng_offdes.design.NbhtD
                NbltD = eng_offdes.design.NbltD

                mbfD  = eng_offdes.design.mbfD
                mblcD = eng_offdes.design.mblcD
                mbhcD = eng_offdes.design.mbhcD
                mbhtD = eng_offdes.design.mbhtD
                mbltD = eng_offdes.design.mbltD

                pifD  = eng_offdes.design.pifD
                pilcD = eng_offdes.design.pilcD
                pihcD = eng_offdes.design.pihcD
                pihtD = eng_offdes.design.pihtD
                piltD = eng_offdes.design.piltD

                if (initializes_engine)
                        #------ force TFOPER to initialize these state variables
                        mbf = 0.0
                        mblc = 0.0
                        mbhc = 0.0
                        pif = 0.0
                        pilc = 0.0
                        pihc = 0.0
                        pt5 = 0.0
                        M2 = 1.0
                        M25 = 1.0
                else
                        #------ use existing state variables as initial guesses
                        mbf = pare[iembf]
                        mblc = pare[iemblc]
                        mbhc = pare[iembhc]
                        pif = max(pare[iepif], 1.1)
                        pilc = max(pare[iepilc], 1.1)
                        pihc = max(pare[iepihc], 1.1)
                        pt5 = pare[iept5]
                        M2 = pare[ieM2]
                        M25 = pare[ieM25]
                end

                Fe = 0.0
                Tt4 = pare[ieTt4]

                if opt_calc_call == CalcMode.FixedTt4OffDes
                        #------ specified Tt4 -- Fe will be computed
                        nothing; #nothing special is done
                elseif opt_calc_call == CalcMode.FixedFeOffDes
                        #------ specified Fe -- Tt4 will be computed (set initial guess here)
                        Fe = pare[ieFe]
                end

                if (Lprint)
                        println(cplab[ip], opt_calc_call, Tt4, Fe)
                        println("Calling TFOPER...", opt_calc_call, opt_cooling, ip)
                        println(DAwsurf, para[iagamV])
                        println(rho0, u0, parg[igWMTO])
                        println(mcore, M2, M25)
                        println("Phiinl, Kinl", Phiinl, Kinl)

                end

                TSFC, Fsp, hfuel, ff,
                Fe, mcore,
                pif, pilc, pihc,
                mbf, mblc, mbhc,
                Nbf, Nblc, Nbhc,
                Tt0, ht0, pt0, cpt0, Rt0,
                Tt18, ht18, pt18, cpt18, Rt18,
                Tt19, ht19, pt19, cpt19, Rt19,
                Tt19c, ht19c, pt19c, cpt19c, Rt19c,
                Tt2, ht2, pt2, cpt2, Rt2,
                Tt21, ht21, pt21, cpt21, Rt21,
                Tt25, ht25, pt25, cpt25, Rt25,
                Tt25c, ht25c, pt25c, cpt25c, Rt25c,
                Tt3, ht3, pt3, cpt3, Rt3,
                Tt4, ht4, pt4, cpt4, Rt4,
                Tt41, ht41, pt41, cpt41, Rt41,
                Tt45, ht45, pt45, cpt45, Rt45,
                Tt49, ht49, pt49, cpt49, Rt49,
                Tt5, ht5, pt5, cpt5, Rt5,
                Tt7, ht7, pt7, cpt7, Rt7,
                u0,
                T2, u2, p2, cp2, R2, M2,
                T25, u25, p25, cp25, R25, M25,
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
                        pifD, pilcD, pihcD, pihtD, piltD,
                        mbfD, mblcD, mbhcD, mbhtD, mbltD,
                        NbfD, NblcD, NbhcD, NbhtD, NbltD,
                        A2, A25, A5, A7,
                        opt_calc_call,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt9, pt9,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Fe,
                        M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
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

                eng_offdes.st18.Tt  = Tt18; eng_offdes.st18.ht  = ht18
                eng_offdes.st18.pt  = pt18; eng_offdes.st18.cpt = cpt18; eng_offdes.st18.Rt = Rt18

                eng_offdes.st19.Tt  = Tt19; eng_offdes.st19.ht  = ht19
                eng_offdes.st19.pt  = pt19; eng_offdes.st19.cpt = cpt19; eng_offdes.st19.Rt = Rt19

                # st19c (PreCoolerOut) — returned by tfoper! but NOT in pare
                eng_offdes.st19c.Tt  = Tt19c; eng_offdes.st19c.ht  = ht19c
                eng_offdes.st19c.pt  = pt19c; eng_offdes.st19c.cpt = cpt19c; eng_offdes.st19c.Rt = Rt19c

                eng_offdes.st2.Tt  = Tt2;  eng_offdes.st2.ht  = ht2
                eng_offdes.st2.pt  = pt2;  eng_offdes.st2.cpt = cpt2;  eng_offdes.st2.Rt = Rt2
                eng_offdes.st2.Ts  = T2;   eng_offdes.st2.ps  = p2
                eng_offdes.st2.Rs  = R2;   eng_offdes.st2.cps = cp2;   eng_offdes.st2.u  = u2
                eng_offdes.st2.A   = A2;   eng_offdes.st2.mdot = mcore

                eng_offdes.st21.Tt  = Tt21; eng_offdes.st21.ht  = ht21
                eng_offdes.st21.pt  = pt21; eng_offdes.st21.cpt = cpt21; eng_offdes.st21.Rt = Rt21

                eng_offdes.st25.Tt  = Tt25; eng_offdes.st25.ht  = ht25
                eng_offdes.st25.pt  = pt25; eng_offdes.st25.cpt = cpt25; eng_offdes.st25.Rt = Rt25
                eng_offdes.st25.Ts  = T25;  eng_offdes.st25.ps  = p25
                eng_offdes.st25.Rs  = R25;  eng_offdes.st25.cps = cp25;  eng_offdes.st25.u  = u25
                eng_offdes.st25.A   = A25

                # st25c (InterCoolerOut) — returned by tfoper! but NOT in pare
                eng_offdes.st25c.Tt  = Tt25c; eng_offdes.st25c.ht  = ht25c
                eng_offdes.st25c.pt  = pt25c; eng_offdes.st25c.cpt = cpt25c; eng_offdes.st25c.Rt = Rt25c

                eng_offdes.st3.Tt  = Tt3;  eng_offdes.st3.ht  = ht3
                eng_offdes.st3.pt  = pt3;  eng_offdes.st3.cpt = cpt3;  eng_offdes.st3.Rt = Rt3

                eng_offdes.st4.Tt  = Tt4;  eng_offdes.st4.ht  = ht4
                eng_offdes.st4.pt  = pt4;  eng_offdes.st4.cpt = cpt4;  eng_offdes.st4.Rt = Rt4

                eng_offdes.st41.Tt  = Tt41; eng_offdes.st41.ht  = ht41
                eng_offdes.st41.pt  = pt41; eng_offdes.st41.cpt = cpt41; eng_offdes.st41.Rt = Rt41

                eng_offdes.st45.Tt  = Tt45; eng_offdes.st45.ht  = ht45
                eng_offdes.st45.pt  = pt45; eng_offdes.st45.cpt = cpt45; eng_offdes.st45.Rt = Rt45

                eng_offdes.st49.Tt  = Tt49; eng_offdes.st49.ht  = ht49
                eng_offdes.st49.pt  = pt49; eng_offdes.st49.cpt = cpt49; eng_offdes.st49.Rt = Rt49

                eng_offdes.st5.Tt  = Tt5;  eng_offdes.st5.ht  = ht5
                eng_offdes.st5.pt  = pt5;  eng_offdes.st5.cpt = cpt5;  eng_offdes.st5.Rt = Rt5
                eng_offdes.st5.Ts  = T5;   eng_offdes.st5.ps  = p5
                eng_offdes.st5.Rs  = R5;   eng_offdes.st5.cps = cp5;   eng_offdes.st5.u  = u5
                eng_offdes.st5.A   = A5

                eng_offdes.st6.Ts  = T6;   eng_offdes.st6.ps  = p6
                eng_offdes.st6.Rs  = R6;   eng_offdes.st6.cps = cp6;   eng_offdes.st6.u  = u6
                eng_offdes.st6.A   = A6

                eng_offdes.st7.Tt  = Tt7;  eng_offdes.st7.ht  = ht7
                eng_offdes.st7.pt  = pt7;  eng_offdes.st7.cpt = cpt7;  eng_offdes.st7.Rt = Rt7
                eng_offdes.st7.Ts  = T7;   eng_offdes.st7.ps  = p7
                eng_offdes.st7.Rs  = R7;   eng_offdes.st7.cps = cp7;   eng_offdes.st7.u  = u7
                eng_offdes.st7.A   = A7

                eng_offdes.st8.Ts  = T8;   eng_offdes.st8.ps  = p8
                eng_offdes.st8.Rs  = R8;   eng_offdes.st8.cps = cp8;   eng_offdes.st8.u  = u8
                eng_offdes.st8.A   = A8

                eng_offdes.st9.u   = u9;   eng_offdes.st9.A   = A9

                # Cooling — mirror the sizing EXIT block (lines 376-378) so that
                # engine_state_to_pare! can write cooling data to pare for both
                # operating modes (FixedCoolingFlowRatio and FixedTmetal).
                # epsrow/Tmrow are populated by tfoper! in place; mofft/mcore
                # are available as tfoper! return values.
                eng_offdes.design.epsrow = SVector{4,Float64}(epsrow[1], epsrow[2], epsrow[3], epsrow[4])
                eng_offdes.design.Tmrow  = SVector{4,Float64}(Tmrow[1],  Tmrow[2],  Tmrow[3],  Tmrow[4])
                eng_offdes.design.fc     = (1.0 - mofft / mcore) * sum(epsrow)

                # Performance rollup scalars (tasopt-j9l.52)
                # BPR is computed from off-design corrected flows (not the design-point input).
                # Fe is the per-engine net thrust (OUTPUT in FixedTt4; INPUT retained in FixedFe).
                eng_offdes.BPR   = mbf / mblc * sqrt(Tt19c / Tt2) * pt2 / pt19c
                eng_offdes.Fe    = Fe
                eng_offdes.TSFC  = TSFC
                eng_offdes.Fsp   = Fsp
                eng_offdes.mfuel = ff * mcore * neng

                # Project thermodynamic state back to pare.  Writes ambient
                # scalars, station fields, and performance rollup scalars.
                # Does NOT write M2/M25/Tt4 (handled below) or design-point
                # scalars (read-only for off-des).
                engine_state_to_pare!(eng_offdes, pare)

                if (Lprint)
                        println("exited TFOPER", Lconv)
                end

                if (!Lconv)
                        #@warn "Convergence failed on operating point: $ip"
                        pare[ieConvFail] = 1.0 #Store convergence failure
                end

                fo = mofft / mcore

                Nf = Nbf * sqrt(Tt2 / Tref)
                N1 = Nblc * sqrt(Tt19c / Tref)
                N2 = Nbhc * sqrt(Tt25c / Tref)

                pare[ieM2] = M2
                pare[ieM25] = M25

                # Fe: mode-dependent (tasopt-j9l.52).
                # FixedTt4OffDes: Fe is a COMPUTED OUTPUT → write bare (engine_state_to_pare!
                #   does NOT write ieFe; see its docstring for the Tt4/Fe analogy).
                # FixedFeOffDes:  Fe is an INPUT (required thrust from mission optimization)
                #   → do NOT overwrite pare[ieFe]; writing Fe_returned (which may differ by
                #   Newton-convergence tolerance) would corrupt the required-thrust value.
                # pare[ieBPR] removed (tasopt-j9l.52): always computed output, written by
                #   engine_state_to_pare! via eng.BPR for both FixedTt4 and FixedFe.
                if opt_calc_call == CalcMode.FixedTt4OffDes
                        pare[ieFe] = Fe
                elseif opt_calc_call == CalcMode.FixedFeOffDes
                        # Tt4 is a COMPUTED OUTPUT in FixedFeOffDes mode; engine_state_to_pare!
                        # does NOT write ieTt4 (Tt4 is a sizing INPUT in the general case).
                        pare[ieTt4] = Tt4
                end

                # println("exited TFOPER call")

        end
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        # Cooling (tasopt-j9l.20):
        # FixedCoolingFlowRatio — Tmrow bare-write removed; now written via
        #   engine_state_to_pare! from eng.design.Tmrow (set in EXIT blocks).
        # FixedCoolingFlowRatio — epsrow is an INPUT (read from pare at entry);
        #   writing it back via engine_state_to_pare! is idempotent.
        # FixedTmetal — epsrow bare-write removed; now written via
        #   engine_state_to_pare! from eng.design.epsrow (set in EXIT blocks).
        # FixedTmetal — fc (iefc) is a COMPUTED OUTPUT only in this mode;
        #   pare[iefc] is NOT touched in FixedCoolingFlowRatio mode (retains
        #   initialization value), so we keep an explicit write here.
        if opt_cooling == CoolingOpt.FixedTmetal
                pare[iefc] = (1.0 - fo) * sum(epsrow)
        end

        # pare[ieTSFC] removed (tasopt-j9l.52): written by engine_state_to_pare! via eng.TSFC
        # pare[ieFsp]  removed (tasopt-j9l.52): written by engine_state_to_pare! via eng.Fsp
        # pare[iemcore] removed (tasopt-j9l.52): written by engine_state_to_pare! via eng.st2.mdot
        pare[iehfuel] = hfuel
        pare[ieff] = ff
        pare[iemofft] = mofft
        pare[iePofft] = Pofft
        pare[iePhiinl] = Phiinl
        pare[ieKinl] = Kinl

        pare[ieNf] = Nbf * sqrt(Tt2 / Tref)
        pare[ieN1] = Nblc * sqrt(Tt19c / Tref)
        pare[ieN2] = Nbhc * sqrt(Tt25c / Tref)
        pare[ieNbf] = Nbf
        pare[ieNblc] = Nblc
        pare[ieNbhc] = Nbhc
        pare[iembf] = mbf
        pare[iemblc] = mblc
        pare[iembhc] = mbhc
        pare[iepif] = pif
        pare[iepilc] = pilc
        pare[iepihc] = pihc

        # Stations and performance scalars written by engine_state_to_pare! (tasopt-j9l.16-52):
        # st0/st18/st19 total state: written via eng.st0/st18/st19
        # st2 total (Tt2/ht2/pt2/cpt2/Rt2): written via eng.st2
        # st21 total (Tt21/ht21/pt21/cpt21/Rt21): written via eng.st21
        # st25 total (Tt25/ht25/pt25/cpt25/Rt25): written via eng.st25
        # st3 total (Tt3/ht3/pt3/cpt3/Rt3): written via eng.st3
        # st4 derived total (ht4/pt4/cpt4/Rt4): written via eng.st4
        #   Tt4 is an INPUT in sizing; written directly above only in FixedFeOffDes mode.
        # st41 total (Tt41/ht41/pt41/cpt41/Rt41): written via eng.st41
        # st45 total (Tt45/ht45/pt45/cpt45/Rt45): written via eng.st45
        # st49 total (Tt49/ht49/pt49/cpt49/Rt49): written via eng.st49
        # st49c is NOT in pare (no ie* indices defined); remains in eng.st49c only.
        # st5 total (Tt5/ht5/pt5/cpt5/Rt5): written via eng.st5
        # st5 static (p5/T5/R5/cp5/u5) and A5: written via eng.st5
        # st6 static (p6/T6/R6/cp6/u6) and A6: written via eng.st6
        # st7 total (Tt7/ht7/pt7/cpt7/Rt7): written via eng.st7
        # st7 static (p7/T7/R7/cp7/u7) and A7: written via eng.st7
        # st8 static (p8/T8/R8/cp8/u8) and A8: written via eng.st8
        # st9 (u9/A9): written via eng.st9
        # st2 static (p2/T2/R2/cp2/u2): written via eng.st2
        # st25 static (p25/T25/R25/cp25/u25): written via eng.st25
        # ieu0: written via eng.st0.u
        # ieTSFC/ieFsp: written via eng.TSFC/eng.Fsp (tasopt-j9l.52)
        # ieFe: conditional bare write above (FixedTt4OffDes only); FixedFeOffDes retains
        #   the required thrust from mission optimization (tasopt-j9l.52)
        # ieBPR: written via eng.BPR (tasopt-j9l.52)
        # iemcore: written via eng.st2.mdot (tasopt-j9l.52)
        # iemfuel: written via eng.mfuel (tasopt-j9l.52)

        pare[ieepf] = epf
        pare[ieeplc] = eplc
        pare[ieephc] = ephc
        pare[ieepht] = epht
        pare[ieeplt] = eplt

        pare[ieetaf] = etaf
        pare[ieetalc] = etalc
        pare[ieetahc] = etahc
        pare[ieetaht] = etaht
        pare[ieetalt] = etalt

        # pare[iemfuel] removed (tasopt-j9l.52): written by engine_state_to_pare! via eng.mfuel

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

        if (Lprint)
                println(" exiting TFCALC")
                println("Tt3 Tt4 u0 u6 u8 fo fc", Tt3, Tt4, u0, u6, u8, fo, pare[iefc])
        end
        
        return ichoke5, ichoke7
end # tfcalc

function check_engine_convergence_failure(pare)
        if sum(pare[ieConvFail, :]) > 0.0 #If any operating point failed to converge
                warn("Some engine points failed to converge: "*string(pare[ieConvFail,:]))
        # else
        #         nothing #All operating points converged
        end
end
