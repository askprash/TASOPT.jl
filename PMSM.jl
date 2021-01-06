"""
PMSM sizes the electric machine for the motor and generator 
    of the turbo electric powertrain

Inputs:
 - P        :   power [W]
 - N        :   Rotational speed [1/s] 
 - ratSplit :   dRot/dStator 
 - σAg      :   Air-gap shear stress


Outputs:
 - Preq   : Required input power
 - η      : Overall efficiency P/Preq
 - W      : Weight of machine [N]
 - l      : Length of machine [m]
 - d      : Diameter of machine [m]


Originally based of code by W. Enders
Modified and updated by P. Prashanth 
"""
function PMSM(P::Float64, ratAsp::Float64, σAg::Float64, ratSplit::Float64, parte::Array{Float64, 1} )

    # Setup/ assumptions
        Vtip = 200 # [m/s] Rotor tip speed with retaining sleeve

        # Rotor dimensions (See Hendershot p89)
            # |T| = |F*r| = (π*D*L)*σ * D/2
            # P = Tω =  (π*D*L)*σ * D/2 *V/(D/2) 
            lRot = sqrt(P/(π*ratAsp*σAg*Vtip))
            dRot = lRot*ratAsp

            N = Vtip/(dRot*π)

        # [TODO] From Wilhelm's aircraft sizing code... but need to check why 240  
        p     = parte[ite_p]

        κM = 240/(2p)
   
    # Rotational speed
        ω = 2π*N    

    # Unpack paramters
        ratAg = parte[ite_ratAg]
        ratM  = parte[ite_ratM ]
        hRS   = parte[ite_hRS  ]
        ratSp = parte[ite_ratSp]
        ratSM = parte[ite_ratSM]
        ratSd = parte[ite_ratSd]
        ratW  = parte[ite_ratW ]
        Nshrt = parte[ite_Nshrt]
        kpf   = parte[ite_kpf  ]

        
        z     = parte[ite_z]
        Br    = parte[ite_Br]

        μ0    = parte[ite_mu0]

        ρMag   = parte[ite_rhoMag  ]
        ρCu    = parte[ite_rhoCu   ]
        ρFe    = parte[ite_rhoFe   ]
        ρSteel = parte[ite_rhoSteel]

        ratShft  = parte[ite_ratShft]
        τMax     = parte[ite_tauMax ]

        ψ  = parte[ite_psi]
    # ------------------
    # Machine geometry
    # ------------------
        Q = P/ω
        Qmax = 1.5*Q

        rRot = dRot/2.0

        hAg = ratAg*dRot + hRS # Air-gap Thickness
        hM  = ratM*hAg

        # Electric frequency
            rRoti = dRot/2 - (hM + hRS) #inner radius of rotor
            f = p * N
            Ω = f * 2π
            Vtip ≤ (rRoti + 0.5*hM)*ω ? println("Warning Vtip test not passed") :
        
        # Calculate Slot geometry
            NS  = ratSp * p  # Number of slots from slots-polepair ratio ratSp
            hS  = ratSM * hM # Slot height
            hSd = ratSd * hS # Slot depression height
            wST = 2π/NS * (rRoti + hM + hAg + hSd + 0.5*hS) # Slot pitch

            wT = wST * ratW # Tooth width
            wS = wST - wT

            κS  = wS/(rRoti + hAg - hRS + hSd + 0.5*hS) # Angular width of stator gaps see W. Enders thesis
            wSi = κS * (rRoti + hAg + hM + hSd)         # Arc length of inner part of stator

            δ = wS/wST

            # Stator back Iron height
                hSBI = rRot/ratSplit - (rRot + hAg + hRS + hSd + hS)

            # Slots/pole/phases
                m   = NS/(2*p*z)  # Acc to Hanselman p125, this should usually be ≤ 2
                NSz = NS/z        # Number of slots per phase

                NSfpc = floor(m*z)
                NSc   = NSfpc - Nshrt  # Actual number of slots per pole
                cSp   = NSc/NS         # Short pitch correction factor

            # End-winding lengths
            # Derived from (Ofori-Tenkorang1996), considering short pitch term 1/2p -> Nsct/Ns
                l_ax  =  π * cSp * (rRoti + hAg + hM + hSd + 0.5*hS) * δ/sqrt(1-δ^2)
                l_cir = 2π * cSp
                
                l_ewind = 2. * sqrt(l_ax^2 + (l_cir/2)^2)
            # Minimum magnet skew angle to prevent Cogging Torque #[TODO] Check this - seems fishy
                κMs = (NSc*2p)/ lcm(Int(NSc*2p), Int(2p)) *1/p * 180.0/π
            
        # Calculate correction factors
            # Calculate pitch factor kp
                α = π * NSc/NSfpc # Angular dispalcement of two sides of the coil
                kp = sin(α/2)
            # Calculation of breadth factor kb
                γ  = 2π * p /NS
                kb = sin(m*γ/2)/(m*sin(γ/2))
            # Winding factor
                kw = kp*kb

            # Skew factor (Jagiela2013)
                κMs_rad = κMs * π/180.
                ks = sin(p*κMs_rad/2) / (p*κMs_rad/2)

            # Magentic gap factor
                rs = rRoti + hM + hAg
                r2 = rRoti + hM
                kg = ((rRoti^(p-1))/(rs^(2*p) - rRoti^(2*p))) * 
                    (  (p/(p+1))            * (r2^(p+1) - rRoti^(p+1)) + 
                        (p/(p-1)) * rs^(2*p) * (rRoti^(1-p) - r2^(1-p))  )
            
        # Calcualte back EMF
            # Accounting for slots, reluctance and leakage
            Kc = 1/(1-(1/((wST/wS)*((5*hAg/wS)+1))));
            hAgEff = Kc*hAg            # Effective Air-gap 
            Cphi   = (p*κM)/180.0      # Flux concentration factor
            K1 = 0.95                  # Leakage factor
            Kr = 1.05                  # Reluctance factor
            μrec = 1.05                # Recoil permeability

            PC  = hM/(hAgEff * Cphi);  # Permeance coefficient
            BAg = ((K1 * Cphi)/(1 + (Kr * μrec/PC) ))*Br;  # Flux density in the air-gap

        # Calculate magnetic flux and internal voltage
            κMrad = κM*(π/180);

            BC1 = (4/π)*BAg*kg*sin(p*κMrad/2);
            λ = 2*rs*lRot*NSz*kw*ks*BC1/p;

            Erms  = Ω*λ/sqrt(2); # RMS back voltage
    
        # Calculation of inductance and total reactance 
            # Air-gap inductance [See 6.685 notes]
                LAg = (z)*(2/π)*(μ0*NSz^2*kw^2*lRot*rs)/(p^2*(hAg+hM));

            # Slot leakage inductance
                Cperm = μ0*((1/3)*(hS/wSi) + hSd/wSi);
                Las = 2*p*lRot*Cperm*(4*(m-Nshrt)+2*Nshrt);
                Lam = 2*p*lRot*Nshrt*Cperm;
                Ls = Las - Lam; # equ. for 3 phases only

            # End-turn inductance 
                areaS = wS*hS; # Slot area
                Le = 0.25*(μ0*wST*NSz^2)*log(wST*sqrt(π)/sqrt(2*areaS));
                #wST only true, if coils are placed in neighboring slots 


            # Total inductance and reactance per phase
                Ltot = LAg + Ls + Le;
                Xtot = Ω*Ltot;


    ## ------------Calculation of machine dimensions and weights---------------
        #Armature
            #Total Armature length per phase (assuming two coils per slot)
                lArm=2*NSz*(lRot + l_ewind); 
            # Armature conductor area
                areaArm = 0.5*areaS*kpf;
            # Total mass of armature conductor #mass=pha*l*area*rho
                mArm = z*lArm*areaArm*ρCu;
        
        
        #Iron /Stator Core
            wSd = parte[ite_wSd]

            # SBI inside radius
                rSBIi = rRoti + hM + hAg + hSd + hS;
            # SBI outside radius
                rSBIo = rSBIi + hSBI;
            # Core mass
                mSBI   = pi*(rSBIo^2 - rSBIi^2)*lRot*ρFe; # SBI
                mTeeth = (NS*wT*hS + 2*pi*(rRoti+hAg)*hSd - NS*hSd*wSd)*lRot*ρFe; # Teeth
                mIron  = mSBI + mTeeth;
        
        # Magnet mass
            mM = (p*κMrad)*((rRoti+hM)^2-rRoti^2)*lRot*ρMag; 
        
        # Shaft mass (Hollow shaft)
            #tauMax=Qmax/Wt Wt=section modulus (Widerstandsmoment)
            #Wt=pi/16*(da^4-di^4)/da, (Dankert,Technische Mechanik, 2018)
            lShft  = 1.2*lRot               #A ssumption to account for bearings, endcaps, etc.
            rShfto = 0.5*((16*Qmax)/(τMax*pi*(1-ratShft^4)))^(1/3) # Outer shaft radius

            mShft  = (π*rShfto^2*(1-ratShft^2))*lShft*ρSteel       # Mass of shaft
            tShft  = rShfto*(1-ratShft)                            # thickness of shaft, just for comparison
        
        
        kServ = parte[ite_kServ]
        # Total mass
            mPMSM = kServ*(mIron + mShft + mM + mArm)
            W = mPMSM * gee
        
        
        # Final machine dimensions
        lPMSM = lRot + 2l_ax # Total length
        dPMSM = 2*rSBIo      # Total diameter without housing
    
    # --------------------------
    # Machine performance
    # --------------------------
        ## Calculate design current and Armature resistance 
            #Armature RMS Current
                #Irms=POutD/(z*cos(ψ)*Erms);
                Irms = 1/(sqrt(2)*3)*Q/(kw*ks*NSz*BAg*(rRoti+hM)*lRot); # (Lipo2017)

            #Armature resistance per phase
                θCu  = parte[ite_thetaCu]
                σCu  = parte[ite_sigCu]
                Tarm = parte[ite_Tarm]

                Rarm = lArm*(1+ θCu*(Tarm - 293.15))/(σCu*areaArm) # Pyrhönen2008 
        
        ##Terminal Voltage
            VaD = sqrt(Erms^2-((Xtot+Rarm)*Irms*cos(ψ))^2)-(Xtot+Rarm)*Irms*sin(ψ)
        
        
        ## Loss Calculations
            #Copper losses
                PLCu = z*Irms^2*Rarm # Total Design Copper Losses
        
            #Iron losses
                BSat = parte[ite_BSat]
                kst  = parte[ite_kst ]
                pb0  = parte[ite_pb0 ]
                Bb0  = parte[ite_Bb0 ]
                fb0  = parte[ite_fb0 ]
                ϵb   = parte[ite_epsb]
                ϵf   = parte[ite_epsf]

                Bt = (wT + wS)/wT*BAg; #Tooth Flux Density (eq. 8.77,Lipo2017IntroToACDesign)
                if Bt > BSat
                    println("ERROR: Bt > Saturation flux density BIronSat")
                end
        
                Bsbi = BAg*rRoti*π/(2*p*kst*hSBI); #BackIron flux density(Hanselman eq 9.7)
                if Bsbi > BSat
                    println("ERROR: Bsbi > Saturation flux density BIronSat")
                end
            
                PLsbi =   mSBI * pb0 * abs(Bsbi/Bb0)^ϵb * abs(f/fb0)^ϵf; # Design SBI Losses
                PLt   = mTeeth * pb0 * abs(  Bt/Bb0)^ϵb * abs(f/fb0)^ϵf; # Design Teeth Losses
                
                PLiron = PLsbi + PLt # Total Design Core Losses
        
            #Windage loss 
                # [TODO] This uses air as the gas, what if we use H2 vented out of the tank... is it stupidly dangerous? Not if you have a sealed enclosure with no O₂?
                Re = N*2π*rRoti*hAg*ρAir/μAir  # Reynold's number in air gap
                Cf = 0.0725*Re^(-0.2)          # Friction coefficient
                PLwind = Cf*pi*ρAir*(2π*N)^3*rRoti^4*lRot  # Design Windage losses

        # Current density
            JarmD = Irms/areaArm;
            if JarmD > 3e7
                println("ERROR: JaD > 3e7 [A/m^2]")
            end
        
        
        ## Design Efficiency and Current Density Calculation
            # Required power and efficiency
            PL   = PLiron + PLCu + PLwind
            Preq = P + PL
            η    = P/Preq
            
            rpm = 60N  # Output rotational speed

return W, Preq, η, rpm, PL, PLiron/PL, PLCu/PL, PLwind/PL

end

function PMSM(P::Float64, N::Float64, parp)

    # λ = 

    # kw = 
    # ks =
    # NSz = 
    # BAg = 
    # rRot =
    # hM   = 
    # lRot =


    # Calculations
        ω = 2π*N
        Q = P/ω

        f = p * N
        Ω = 2π*f

        Erms = Ω * λ / sqrt(2)
        Irms = 1/(3*sqrt(2)) * Q / 
               (kw * ks * NSz * BAg * (rRot + hM)*lRot)

        Jarm = Irms/areaArm

        # Losses
            PLsbi = mSBI  * pb0 * abs(Bsbi/Bb0)^epsb * abs(f/fb0)^ϵf
            PLt   = mTeeth* pb0 * abs(  Bt/Bb0)^epsb * abs(f/fb0)^ϵf
            PLiron = PLsbi + PLt

            Re = ω*rRoti*hAg*ρair/μair
            Cf = 0.0725/Re^0.2
            PLwind = ρair * π * Cf * ω^3 * rRoti^4 * lRot

            PLcu = z * Irms^2 * Rarm

            PL = PLiron + PLwind + PLcu

        Preq = P + PL
        η    = P/Preq

    return Preq, η, PL
end
