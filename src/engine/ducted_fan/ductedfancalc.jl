"""
    ductedfancalc!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)

Calls function ductedfansize! or ductedfanoper! for one operating point.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `case::String`: case identifier, e.g. "sizing" or "off_design"
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `initializes_engine::Bool`: flag to initialize engine
      - `true`: initialize variables for iteration in engine
      - `false`: use current variables as initial guesses in engine
    - `iterw::Int64`: sizing loop iteration

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function ductedfancalc!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)
    #Unpack data storage arrays
    parg, _, para, pare, options, _, _, wing, _, _, eng, _ = unpack_ac(ac, imission, ip=ip)
    iBLIc = eng.model.has_BLI_cores

    neng = parg[igneng]
    S = wing.layout.S

    # Typed per-point engine state for this mission/point (tasopt-j9l.45.3)
    eng_ip = ac.missions[imission].points[ip].engine

    Fsp   = eng_ip.Fsp
    pif   = eng_ip.pif
    pid   = eng_ip.design.pid
    pifn  = eng_ip.design.pifn

    epolf = eng_ip.design.epolf

    M2  = eng_ip.design.M2
    M0  = eng_ip.M0
    Tt0 = eng_ip.st0.Tt
    ht0 = eng_ip.st0.ht
    pt0 = eng_ip.st0.pt
    cpt0 = eng_ip.st0.cpt
    Rt0 = eng_ip.st0.Rt
    p0   = eng_ip.p0
    a0   = eng_ip.a0
    rho0 = eng_ip.rho0
    mu0  = eng_ip.mu0
    T0   = eng_ip.T0
    u0   = eng_ip.st0.u

    #Heat exchanger variables
    Δh_radiator = eng_ip.RadiatorDeltah
    Δp_radiator = eng_ip.RadiatorDeltap

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
            CDAwing = para[iaCDwing] * S
            DAwsurf = CDAwing * (1.0 - fDwake)
            KAwTE = DAwsurf

            fBLIw = parg[igfBLIw]

            #----- set ingested  PKinl-PVinl = Phiinl  for one engine
            Phiinl = 0.5 * rho0 * u0^3 * (DAfsurf * fBLIf + DAwsurf * fBLIw) / neng
            Kinl = 0.5 * rho0 * u0^3 * (KAfTE * fBLIf + KAwTE * fBLIw) / neng
    end
    #- - - - - - - - - - - - - - - - - - - - - - -

    # #--------------------------------------------------------------------------
    if (case == "design")
        #----- engine sizing case

        Fe = eng_ip.Fe  #ducted fan sized for a given thrust

        TSEC, Fsp, Pfan, mfan,
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
        etaf, Lconv = ductedfansize!(gee, M0, T0, p0, a0, M2,
                        Fe, Phiinl, Kinl, iBLIc,
                        pif,
                        pid, pifn,
                        epolf,
                        Δh_radiator, Δp_radiator
                        )

        mbf = mfan * sqrt(Tt2 / Tref) / (pt2 / pref)
        M7 = u7 / sqrt(T7 * R7 * cp7 / (cp7 - R7))

        pifD = pif
        mbfD = mbf

        Nf = 1.0 #Arbitrarily set to 1 as only ratios matter
        Nbf = Nf / sqrt(Tt2 / Tref)
        NbfD = Nbf

        #----- store design-point parameters (dual-write: bare pare + typed state)
        pare[ieA2] = A2;  eng_ip.design.A2  = A2
        pare[ieA7] = A7;  eng_ip.design.A7  = A7

        pare[iembfD] = mbfD;  eng_ip.design.mbfD = mbfD
        pare[iepifD] = pifD;  eng_ip.design.pifD = pifD
        pare[ieNbfD] = Nbf;   eng_ip.design.NbfD = Nbf

    else
        #----- fixed parameters
        A2   = eng_ip.design.A2
        A7   = eng_ip.design.A7

        mbfD = eng_ip.design.mbfD
        pifD = eng_ip.design.pifD
        NbfD = eng_ip.design.NbfD

        if initializes_engine
                #------ initialize these state variables
                mbf = 0.0
                pif = 0.0
                M2 = 1.0
        else
                #------ use existing state variables as initial guesses
                mbf = eng_ip.mbf
                pif = eng_ip.pif
                M2  = eng_ip.design.M2
        end
        if ip in range(ipstatic, iprotate)
            #Power is specified, thrust will be computed -- in takeoff
            Feng = 0.0
            Peng = eng_ip.Pfanmax
            iPspec = true

            TSEC, Fsp, Feng, Pfan, mfan,
            pif, mbf, Nbf,
            Tt0, ht0, pt0, cpt0, Rt0,
            Tt18, ht18, pt18, cpt18, Rt18,
            Tt2, ht2, pt2, cpt2, Rt2,
            Tt21, ht21, pt21, cpt21, Rt21,
            Tt7, ht7, pt7, cpt7, Rt7,
            u0, T2, u2, p2, cp2, R2, M2,
            T7, u7, p7, cp7, R7, M7,
            T8, u8, p8, cp8, R8, M8, A8,
            epf, etaf = ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                            Phiinl, Kinl, iBLIc,
                            pid, pifn,
                            pifD,
                            mbfD, NbfD,
                            A2, A7,
                            epolf,
                            Feng, Peng,
                            M2, pif, mbf,
                            Δh_radiator, Δp_radiator,
                            iPspec)
            pare[ieFe] = Feng;  eng_ip.Fe = Feng

        else #Thrust is specified, power to be computed
            Feng = eng_ip.Fe
            Peng = 0.0
            iPspec = false

            TSEC, Fsp, Feng, Pfan, mfan,
            pif, mbf, Nbf,
            Tt0, ht0, pt0, cpt0, Rt0,
            Tt18, ht18, pt18, cpt18, Rt18,
            Tt2, ht2, pt2, cpt2, Rt2,
            Tt21, ht21, pt21, cpt21, Rt21,
            Tt7, ht7, pt7, cpt7, Rt7,
            u0, T2, u2, p2, cp2, R2, M2,
            T7, u7, p7, cp7, R7, M7,
            T8, u8, p8, cp8, R8, M8, A8,
            epf, etaf = ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                            Phiinl, Kinl, iBLIc,
                            pid, pifn,
                            pifD,
                            mbfD, NbfD,
                            A2, A7,
                            epolf,
                            Feng, Peng,
                            M2, pif, mbf,
                            Δh_radiator, Δp_radiator,
                            iPspec)
        end
    end
    HTRf = parg[igHTRf]
    dfan = sqrt(4.0 * A2 / (pi * (1.0 - HTRf^2)))
    parg[igdfan] = dfan

    pare[iePfan] = Pfan;   eng_ip.Pfan  = Pfan
    pare[ieTSEC] = TSEC;   eng_ip.TSEC  = TSEC
    pare[ieTSFC] = TSEC / 120e6;  eng_ip.TSFC = TSEC / 120e6  #TODO change when fuel is not burnt
    pare[ieFsp]  = Fsp;    eng_ip.Fsp   = Fsp
    pare[iemfan] = mfan;   eng_ip.mfan  = mfan
    pare[iePhiinl] = Phiinl;  eng_ip.Phiinl = Phiinl
    pare[ieKinl]   = Kinl;    eng_ip.Kinl   = Kinl

    pare[iembf] = mbf;  eng_ip.mbf = mbf
    pare[iepif] = pif;  eng_ip.pif = pif
    pare[ieNbf] = Nbf;  eng_ip.Nbf = Nbf
    Nf_val = Nbf * sqrt(Tt2 / Tref)
    pare[ieNf] = Nf_val;  eng_ip.Nf = Nf_val

    pare[ieTt0]  = Tt0;   eng_ip.st0.Tt  = Tt0
    pare[ieht0]  = ht0;   eng_ip.st0.ht  = ht0
    pare[iept0]  = pt0;   eng_ip.st0.pt  = pt0
    pare[iecpt0] = cpt0;  eng_ip.st0.cpt = cpt0
    pare[ieRt0]  = Rt0;   eng_ip.st0.Rt  = Rt0

    pare[ieTt18]  = Tt18;   eng_ip.st18.Tt  = Tt18
    pare[ieht18]  = ht18;   eng_ip.st18.ht  = ht18
    pare[iept18]  = pt18;   eng_ip.st18.pt  = pt18
    pare[iecpt18] = cpt18;  eng_ip.st18.cpt = cpt18
    pare[ieRt18]  = Rt18;   eng_ip.st18.Rt  = Rt18

    pare[ieTt2]  = Tt2;   eng_ip.st2.Tt  = Tt2
    pare[ieht2]  = ht2;   eng_ip.st2.ht  = ht2
    pare[iept2]  = pt2;   eng_ip.st2.pt  = pt2
    pare[iecpt2] = cpt2;  eng_ip.st2.cpt = cpt2
    pare[ieRt2]  = Rt2;   eng_ip.st2.Rt  = Rt2

    pare[ieTt21]  = Tt21;   eng_ip.st21.Tt  = Tt21
    pare[ieht21]  = ht21;   eng_ip.st21.ht  = ht21
    pare[iept21]  = pt21;   eng_ip.st21.pt  = pt21
    pare[iecpt21] = cpt21;  eng_ip.st21.cpt = cpt21
    pare[ieRt21]  = Rt21;   eng_ip.st21.Rt  = Rt21

    pare[ieTt7]  = Tt7;   eng_ip.st7.Tt  = Tt7
    pare[ieht7]  = ht7;   eng_ip.st7.ht  = ht7
    pare[iept7]  = pt7;   eng_ip.st7.pt  = pt7
    pare[iecpt7] = cpt7;  eng_ip.st7.cpt = cpt7
    pare[ieRt7]  = Rt7;   eng_ip.st7.Rt  = Rt7

    pare[ieu0] = u0;  eng_ip.st0.u = u0

    pare[iep2]  = p2;   eng_ip.st2.ps  = p2
    pare[ieT2]  = T2;   eng_ip.st2.Ts  = T2
    pare[ieR2]  = R2;   eng_ip.st2.Rs  = R2
    pare[iecp2] = cp2;  eng_ip.st2.cps = cp2
    pare[ieu2]  = u2;   eng_ip.st2.u   = u2

    pare[iep7]  = p7;   eng_ip.st7.ps  = p7
    pare[ieT7]  = T7;   eng_ip.st7.Ts  = T7
    pare[ieR7]  = R7;   eng_ip.st7.Rs  = R7
    pare[iecp7] = cp7;  eng_ip.st7.cps = cp7
    pare[ieu7]  = u7;   eng_ip.st7.u   = u7

    pare[iep8]  = p8;   eng_ip.st8.ps  = p8
    pare[ieT8]  = T8;   eng_ip.st8.Ts  = T8
    pare[ieR8]  = R8;   eng_ip.st8.Rs  = R8
    pare[iecp8] = cp8;  eng_ip.st8.cps = cp8
    pare[ieu8]  = u8;   eng_ip.st8.u   = u8

    pare[ieA8] = A8;  eng_ip.st8.A = A8

    pare[ieepf]  = epf;   eng_ip.epf  = epf

    pare[ieetaf] = etaf;  eng_ip.etaf = etaf


    if (M7 <= 0.999999)
            ichoke7 = 0
    else
            ichoke7 = 1
    end

end
