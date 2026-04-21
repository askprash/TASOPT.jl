"""
    fly_mission!(ac, imission, itermax, initializes_engine, opt_prescribed_cruise_parameter, printTO)

Runs the aircraft through the specified mission, computing and converging the fuel weight. Formerly, `fly_offdesign_mission!()`.

!!! details "🔃 Inputs and Outputs"
**Inputs:**
- `ac::aircraft`: Aircraft with first mission being the design mission
- `imission::Int64`: Off design mission to run (Default: 1)
- `itermax::Int64`: Maximum iterations for sizing loop
- `initializes_engine::Boolean`: Use design case as initial guess for engine state if true
- `opt_prescribed_cruise_parameter::String`: option for whether cruise altitude or lift coefficient is specified. Options are "altitude" or "lift_coefficient"
- `printTO::Bool`: print takeoff convergence table if true
**Outputs:**
- No explicit outputs. Computed quantities are saved to `par` arrays of `aircraft` model for the mission selected

"""
function fly_mission!(ac, imission = 1; itermax = 35, initializes_engine = true, 
        opt_prescribed_cruise_parameter = "altitude", printTO = true)
    if ~ac.is_sized[1]
        error("Aircraft not sized. Please size the aircraft before running the mission.")
    end
    
    #Extract aircraft components and storage arrays
    parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, eng = unpack_ac(ac, imission)
    
    parad = ac.parad
    pared = ac.pared

    time_propsys = 0.0

    tolerW = 1.0e-8
    errw   = 1.0

    #Warn user if CL or altitude are being unexpectedly overwritten
    if compare_strings(opt_prescribed_cruise_parameter, "altitude")
        if ac.para[iaCL, ipcruise1,1] != ac.para[iaCL, ipcruise1, imission] #if CLs are prescribed but overwritten 
            @warn "An off-design CL is specified, but the analysis is set to rewrite it. Please review your inputs, especially `opt_prescribed_cruise_parameter`, if this is not the desired behavior" maxlog=1   
        end

    elseif compare_strings(opt_prescribed_cruise_parameter, "CL")
        if ac.para[iaalt, ipcruise1,1] != ac.para[iaalt, ipcruise1, imission] #if CLs alts prescribed but overwritten
            @warn "An off-design altitude is specified, but the analysis is set to rewrite it. Please review your inputs, especially `opt_prescribed_cruise_parameter`, if this is not the desired behavior" maxlog=1   
        end
    end
    
#------ mission-varying excrescence factors disabled in this version
#-      ( also commented out in getparm.f )
#        para(iafexcdw,ip) = parm[imfexcdw]
#        para(iafexcdt,ip) = parm[imfexcdt]
#        para(iafexcdf,ip) = parm[imfexcdf]

    #Initialize arrays with the design mission values if desired
    if (initializes_engine)
        #----- use design case as initial guess for engine state
        pare[:,:] .= pared[:,:]  # bare-pare copy for downstream callers not yet migrated to typed state
        # Initialize typed engine state from design mission typed state.
        # Design typed state is authoritative after .14.3 (tasopt-j9l.45.14.3); deepcopy
        # replaces the old bare-pare-copy + pare_to_engine_state! pattern (.14.5).
        des_pts = ac.missions[1].points
        mis_pts = ac.missions[imission].points
        for ip in eachindex(mis_pts)
            mis_pts[ip].engine = deepcopy(des_pts[ip].engine)
        end
    else
        ac.missions[imission].points[ipcruise1].engine.u0 = ac.missions[1].points[ipcruise1].engine.u0
    end

    for ip = ipstatic: ipdescentn
        para[iaCfnace,ip] = parad[iaCfnace,ip]
    end

    #Calculate sea level temperature corresponding to TO conditions
    altTO = parm[imaltTO] 
    T_std = atmos(altTO).T
    ΔTatmos = parm[imT0TO] - T_std #temperature difference such that T(altTO) = T0TO
    parm[imDeltaTatm] = ΔTatmos

    # Calculates surface velocities, boundary layer, wake 
    fuselage_drag!(fuse, parm, para, ipcruise1)
    broadcast_fuselage_drag!(para, ipcruise1) #Broadcast fuselage drag to all flight points

# ===================================================================
# ---- max range and this mission range
    Rangemax = parg[igRange]
    Rangetot = parm[imRange]

#---- max TO weight
    WMTO = parg[igWMTO]

# ---- zero-fuel weight for this mission
    Wzero = WMTO-
          parg[igWfuel]-
          parg[igWpay]+
          parm[imWpay]

# ===================================================================
# ---- initial fuel and gross takeoff weight estimates from Breguet, R ~ ln(1+f)
    gmax = log(1.0 + parg[igWfuel]/Wzero)
    gmaxp = gmax * Rangetot/Rangemax
    Wfuel = (exp(gmaxp) - 1.0) * Wzero
    WTO = Wzero + Wfuel

    parm[imWfuel] = Wfuel
    parm[imWTO]   = WTO

#---- scale initial weight fractions by takeoff and descent weight ratios
    rTO = WTO/WMTO
    rDE = Wzero/(WMTO-parg[igWfuel])

    para[iafracW, ipstatic ] = parad[iafracW,ipstatic ]*rTO
    para[iafracW, iprotate ] = parad[iafracW,iprotate ]*rTO
    para[iafracW, iptakeoff] = parad[iafracW,iptakeoff]*rTO
    para[iafracW, ipcutback] = parad[iafracW,ipcutback]*rTO

    # Climb
    @inbounds for ip = ipclimb1:ipclimbn
          para[iafracW,ip] = parad[iafracW,ip] * rTO
    end
    # Cruise
    @inbounds for ip = ipcruise1:ipcruisen
          frac = float(ip       -ipcruise1)/
                float(ipcruisen-ipcruise1)
          rCR = rTO*(1.0-frac) + rDE*frac
          para[iafracW,ip] = parad[iafracW,ip] * rCR
    end
    # Descent
    para[iafracW,ipdescent1:ipdescentn] .= parad[iafracW,ipdescent1:ipdescentn] .* rDE
    para[iagamV,:] .= parad[iagamV,:]

#---- estimate takeoff speed and set V,Re over climb and descent
#-    (needed to start trajectory integration)
    ip = iptakeoff
    eng_des_to = ac.missions[1].points[ip].engine
    eng_mis_to = ac.missions[imission].points[ip].engine
    VTO = eng_des_to.u0 * sqrt(eng_des_to.rho0 / eng_mis_to.rho0)
    ReTO = VTO * eng_mis_to.rho0 / eng_mis_to.mu0

    ip = ipcruise1
    VCR = ac.missions[1].points[ip].engine.u0
    ReCR = parad[iaReunit,ip]

    for ip = iprotate: ipclimb1
      ac.missions[imission].points[ip].engine.u0 = VTO
      para[iaReunit,ip] = ReTO
    end
    for ip = ipclimb1+1 : ipclimbn
      frac = float(ip-ipclimb1) / float(ipclimbn-ipclimb1)
      V  =  VTO*(1.0-frac) +  VCR*frac
      Re = ReTO*(1.0-frac) + ReCR*frac
      ac.missions[imission].points[ip].engine.u0 = V
      para[iaReunit,ip] = Re
    end
    for ip = ipdescent1: ipdescentn
      frac = float(ip-ipdescent1) / float(ipdescentn-ipdescent1)
      V  =  VTO*frac +  VCR*(1.0-frac)
      Re = ReTO*frac + ReCR*(1.0-frac)
      ac.missions[imission].points[ip].engine.u0 = V
      para[iaReunit,ip] = Re
    end

#--------------------------------------------------------------------------
#---- set wing pitching moment constants
    b  = wing.layout.span
    bs = wing.layout.break_span
    bo = wing.layout.root_span
    sweep = wing.layout.sweep
    Xaxis = wing.layout.spar_box_x_c
    λs = wing.inboard.λ
    λt = wing.outboard.λ
    AR = wing.layout.AR
    fLo =  wing.fuse_lift_carryover
    fLt =  wing.tip_lift_loss

    ip = iptakeoff
    cmpo = para[iacmpo,ip]
    cmps = para[iacmps,ip]
    cmpt = para[iacmpt,ip]
    γt = wing.outboard.λ*para[iarclt,ip]
    γs = wing.inboard.λ*para[iarcls,ip]

    CMw0, CMw1 = wing_CM(b, bs, bo, sweep, Xaxis,
                            λt,λs,γt,γs, 
                            AR,fLo,fLt,cmpo,cmps,cmpt)

    para[iaCMw0, ipstatic:ipclimb1] .= CMw0
    para[iaCMw1, ipstatic:ipclimb1] .= CMw1

    ip = ipcruise1
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]

    γt = wing.outboard.λ*para[iarclt, ip]
    γs = wing.inboard.λ*para[iarcls, ip]
    
    CMw0, CMw1 = wing_CM(b, bs, bo, sweep, Xaxis,
                      λt,λs,γt,γs, 
                      AR,fLo,fLt,cmpo,cmps,cmpt)
    
    para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
    para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1
    
    ip = ipdescentn
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
    γt = wing.outboard.λ*para[iarclt, ip]
    γs = wing.inboard.λ*para[iarcls, ip]

    CMw0, CMw1 = wing_CM(b, bs, bo, sweep, Xaxis,
                      λt,λs,γt,γs, 
                      AR,fLo,fLt,cmpo,cmps,cmpt)

    para[iaCMw0, ipdescentn] = CMw0
    para[iaCMw1, ipdescentn] = CMw1

#---- tail pitching moment constants
    bh      = htail.layout.span
    boh     = htail.layout.root_span
    sweeph  = htail.layout.sweep
    λh      = htail.outboard.λ
    ARh     = htail.layout.AR
    fLoh = 0.
    fLth = fLt
    cmph = 0.

    CMh0, CMh1 = wing_CM(bh, boh, boh, sweeph, Xaxis, λh, 1.0, λh, 1.0,
    ARh, fLoh, fLth, 0.0, 0.0, 0.0)

    para[iaCMh0, :] .= CMh0
    para[iaCMh1, :] .= CMh1

    # Initialize previous weight iterations
    WTO1, WTO2, WTO3 = zeros(Float64, 3) #1st-previous to 3rd previous iteration weight for convergence criterion

    resetHXs(pare, ac.missions[imission].points) #Reset heat exchanger parameters

#---- no convergence yet
    Lconv = false

# -------------------------------------------------------    
#                   Weight loop
# -------------------------------------------------------    
  @inbounds for  iterw = 1:itermax

    if iterw == itermax
        println("Reached max iterations in weight sizing loop!")
    end

    rlx = 1.0
    if (iterw > itermax-5)
          rlx = 0.5
    end

    # Calculate start-of-cruise altitude or CL from each other by ensuring L=W
    calculate_cruise_altitude_or_CL!(opt_prescribed_cruise_parameter, WMTO, ac, imission)
    
    if !(options.has_wing_fuel) #If fuel is stored in the fuselage
        #Analyze pressure evolution in tank and store the vented mass flow rate
        _, _, _, _, _, _, _, Mvents, _, _ = CryoTank.analyze_TASOPT_tank(ac, fuse_tank.t_hold_orig, fuse_tank.t_hold_dest, imission)
        parm[imWfvent] = Mvents[end] * gee #Store vented weight
    end

    # Calling mission
    time_propsys += _mission_iteration!(ac, imission, false, calculate_cruise = true) #Calculate start of cruise too
    # println(parm[imWfuel,:])

    #Simulate heat exchanger performance if the engine contains any
    if eng.model.model_name == "ducted_fan"
        # dual-write to typed state (tasopt-keh)
        for ip in eachindex(ac.missions[imission].points)
            pt = ac.missions[imission].points[ip]
            pt.engine.RadCoolantT = eng.data.FC_temperature[ip, imission]
            pt.engine.RadCoolantP = eng.data.FC_pressure[ip, imission]
            pt.engine.Qradiator   = eng.data.FC_heat[ip, imission]
        end

    end     
    HXOffDesign!(eng.heat_exchangers, pare, ac.options.ifuel, imission, ac.missions[imission].points)

#-------------------------------------------------------------------------

# Convergence tests
    
    WTO = parm[imWTO]
    errw1 = (WTO - WTO1)/WTO
    errw2 = (WTO - WTO2)/WTO
    errw3 = (WTO - WTO3)/WTO

    errw = max(abs(errw1), abs(errw2), abs(errw3))

    if (errw <= tolerW) 
          Lconv = true

          break
    end

#-----------------------------------------------------------------
#---- set previous-iteration weights for next iteration
    WTO3 = WTO2
    WTO2 = WTO1
    WTO1 = parm[imWTO]

    end

    #Check if all engine points have converged, warns if not
    check_engine_convergence_failure(ac, imission)

    #run takeoff calculation if converged (checks feasibility; populates iprotate entries)
    if Lconv
        takeoff!(ac, imission=imission, printTO=printTO)
    end

return 
end

"""
    calculate_cruise_altitude_or_CL!(opt_prescribed_cruise_parameter, WMTO, ac, imission)

Calculates the cruise altitude or lift coefficient based on the specified option. If "altitude" is selected, it calculates the lift coefficient
from the weight and density. If "CL" is selected, it calculates the altitude from the lift coefficient and updates the fuselage drag.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `opt_prescribed_cruise_parameter::String`: option for whether cruise altitude or lift coefficient is specified. Options are "altitude" or "CL"
    - `WMTO::Float64`: Maximum takeoff weight (N)
    - `ac`: aircraft model object
    - `imission::Int64`: Mission index
**Outputs:**
    - No explicit outputs. Computed quantities are saved to `par` arrays of `aircraft` model for the mission selected
"""
function calculate_cruise_altitude_or_CL!(opt_prescribed_cruise_parameter, WMTO, ac, imission)
    parg, parm, para, _, _, fuse, _, wing, _, _, _ = unpack_ac(ac, imission)

    #Calculate ΔT for the atmosphere
    altTO = parm[imaltTO] 
    T_std = atmos(altTO).T
    ΔTatmos = parm[imT0TO] - T_std #temperature difference such that T(altTO) = T0TO

    eng_cr = ac.missions[imission].points[ipcruise1].engine
    ρ0 = eng_cr.rho0
    ρcab = max(parg[igpcabin], eng_cr.p0) / (RSL * TSL)
    WbuoyCR = (ρcab - ρ0) * gee * parg[igcabVol]

    ip = ipcruise1
    We = WMTO * para[iafracW, ip]
    u0 = eng_cr.u0
    BW = We + WbuoyCR # Weight including buoyancy
    S = wing.layout.S

    if compare_strings(opt_prescribed_cruise_parameter, "altitude")
        CL = BW / (0.5*ρ0*u0^2*S) #Find CL from L=W
        para[iaCL, ipclimb1+1:ipdescentn-1] .= CL #Store CL in climb, cruise and descent phases
    
    elseif compare_strings(opt_prescribed_cruise_parameter, "CL")
        CL = para[iaCL, ip]
        ρ0 = BW / (0.5*u0^2*S*CL) #Find density from L=W
        para[iaalt, ip] = find_altitude_from_density(ρ0, ΔTatmos) #Store altitude

        set_ambient_conditions!(ac, ipcruise1, im = imission)
        #Update fuselage drag for the new altitude
        fuselage_drag!(fuse, parm, para, ipcruise1)
        broadcast_fuselage_drag!(para, ipcruise1) 
    end
end
