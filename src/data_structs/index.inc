#---- indices for missions to be examined 
#-      (different Ranges Payloads ...)
      const imwOpt = 1 
      const imRange = 2
      const imWpay = 3 
      const imWperpax = 4
      const imaltTO = 5
      const imT0TO = 6
      const imWTO = 7
      const imWfuel = 8
      const imPFEI = 9
      const imV1 = 10
      const imV2 = 11
      const imtTO = 12
      const imFTO = 13
      const iml1 = 14
      const imlTO = 15
      const imlBF = 16
      const imlCB = 17
      const imgamVTO = 18
      const imgamVBF = 19
      const imgamVCB = 20
      const imgamVDE1 = 21
      const imgamVDEn = 22
      const imthCB = 23
      const imxCB = 24
      const imzCB = 25
      const imxFO = 26
      const imzFO = 27
      const imdBSL = 28
      const imdBCB = 29
      const imdBFO = 30
      const imfexcdw = 31
      const imfexcdt = 32
      const imfexcdf = 33
      const imDeltaTatm = 34
      const imWfvent = 35
      const imtotal = 35

#---- indices for geometry (airframe) variables and other sizing variables
      const igRange    =   1
      const igWMTO     =   2
      const igWpay     =   3
      const igWfuel    =   4
      const igWfmax    =   5
      const igrWfmax   =   6
      const igWebare   =  7
      const igWnace    =  8
      const igWeng     =  9
      const igdxWfuel  =  10
      const igxCGfwd   =  11
      const igxCGaft   =  12
      const igfreserve =  13
      const igfeadd    =  14
      const igfpylon   =  15
      const igNlift    =  16
      const igNland    =  17
      const igVne      =  18
      const igneng     =  19
      const igGearf    =  20
      const igfTt4CL1  =  21
      const igfTt4CLn  =  22
      const igHTRf     =  23
      const igHTRlc    =  24
      const igHTRhc    =  25
      const igrSnace   =  26
      const igrVnace   =  27
      const igfSnace   =  28
      const igpcabin   =  29
      const igdeltap   =  30
      const igxeng     =  31
      const igyeng     = 32
      const igcabVol   = 33
      const igrpayfwd  = 34
      const igrpayaft  = 35
      const igxNP      = 36
      const igCMVf1    = 37
      const igCLMf0    = 38
      const igdCLnda   = 39
      const igdCLhdCL  = 40
      const igdCLndCL  = 41
      const igCLhspec  = 42
      const igCLveout  = 43
      const igfCDhcen  = 44
      const igrMh      = 45
      const igrMv      = 46
      const igsigfac   = 47
      const igrhofuel  = 48
      const igfBLIw    = 49
      const igfBLIf    = 50
      const igdfan     = 51
      const igdlcomp   = 52
      const igdhcomp   = 53
      const iglnace    = 54
      const igA5       = 55
      const igA7       = 56
      const igTmetal   = 57
      const igcdefan   = 58
      const igCDgear   = 59
      const igCDspoil  = 60
      const igmuroll   = 61
      const igmubrake  = 62
      const ighobst    = 63
      const iglBFmax   = 64
      const iggtocmin  = 65
      const igdBSLmax  = 66
      const igdBCBmax  = 67
      const igmofWpay  = 68
      const igmofWMTO  = 69
      const igPofWpay  = 70
      const igPofWMTO  = 71

      # Adding indices to store TE weights and locations
      const igWtshaft  = 72
      const igWgen     = 73
      const igWinv     = 74
      const igWmot     = 75
      const igWfan     = 76
      const igWftank   = 77

      const igxtshaft  = 78
      const igxgen     = 79
      const igxinv     = 80
      const igxmot     = 81
      const igxfan     = 82
      const igxftank   = 83
      
      const igxcables  = 84
      const igWcables  = 85

      const igxcat     = 86
      const igWcat     = 87

      const igWtesys   = 88
      const igxWtesys  = 89

      const iglftank   = 90
      const igWinsftank= 91
      const igxWftank  = 92
      
      const igRftank   = 93
      const igWc3des   = 94

      const igdaftfan  = 95
      const iglnaceaft = 96

      const igneout    = 97
      const igyeout    = 98
      const igyeinn    = 99
      const iglftankin = 100 #lftank input
      const igLHVfuel  = 101
      const igWfburn   = 102

      const igWaftfan  = 103
      const igWfanGB   = 104
      const igWaftfanGB   = 105
      const igWrect    = 106
      const igWtms     = 107

      const igfanPCT   = 108

      const iggamVtarget = 109

      const igWpaymax    = 110

      const ig_YC = 111
      const ig_YH = 112
      const ig_YO = 113
      const ig_YN = 114
      const ig_MW = 115
      const ighfuel      = 116
      const igTfuel      = 117

      const igxftankaft  = 118
      const igdxeng2wbox = 119

      const igxWfuel = 120

      const igWfvent = 121
      const igWHXs = 122
      const igfeng = 123
      const igtotal    = 123

#---- indices for aero variables at each mission point
      const iaalt     =  1
      const iatime    =  2
      const iaRange   =  3
      const iafracW   =  4
      const iaWbuoy   =  5
      const iaMach    =  6
      const iaReunit  =  7
      const iagamV    =  9
      const iaCL      = 10
      const iaCD      = 11
      const iaCDi     = 12
      const iaDAfsurf = 13
      const iaDAfwake = 14
      const iaKAfTE   = 15
      const iaPAfinf  = 16
      const iaCDfuse  = 17
      const iaCDwing  = 18
      const iaCDhtail = 19
      const iaCDvtail = 20
      const iaCDnace  = 21
      const iaCDstrut = 22
      const iaCDover  = 23
      const iaCfnace  = 24
      const iarcls    = 25
      const iarclt    = 26
      const iaclpo    = 27
      const iaclps    = 28
      const iaclpt    = 29
      const iacmpo    = 30
      const iacmps    = 31
      const iacmpt    = 32
      const iafduo    = 33
      const iafdus    = 34
      const iafdut    = 35
      const iaclpmax  = 36
      const iaCLh     = 37
      const iaCMw0    = 38
      const iaCMw1    = 39
      const iaCMh0    = 40
      const iaCMh1    = 41
      const iaCMwing  = 42
      const iaCMfuse  = 43
      const iaCMtail  = 44
      const iaxCG     = 45
      const iaxCP     = 46
      const iaxNP     = 47
      const iafexcdw  = 48
      const iafexcdt  = 49
      const iafexcdf  = 50
      const iacdfw    = 51
      const iacdpw    = 52
      const iaspaneff = 53
      const iaRerefw  = 54
      const iaaRexp   = 55
      const iacdft    = 56
      const iacdpt    = 57
      const iaRereft  = 58
      const iacdfs    = 59
      const iacdps    = 60
      const iaRerefs  = 61
      const iadCDBLIf = 62
      const iadCDBLIw = 63
      const iaROC     = 64
      const iaROCdes  = 65
      const iatotal   = 65

#---- max number of HPT blade rows (epsrow/Tmrow vector length)
ncrowx = 4

#---- indices for indexing mission points
      const ipstatic = 1
      const iprotate = 2
      const iptakeoff = 3
      const ipcutback = 4
      const ipclimb1 = 5
      const ipclimb2 = 6
      const ipclimb3 = 7
      const ipclimb4 = 8
      const ipclimb5 = 9
      const ipcruise1 = 10
      const ipcruise2 = 11
      const ipdescent1 = 12
      const ipdescent2 = 13
      const ipdescent3 = 14
      const ipdescent4 = 15
      const ipdescent5 = 16
      const iptest = 17

      # generate label array programmatically from ip vars in TASOPT module namespace
      # if it doesn't already exist
      if !isdefined(@__MODULE__, :ip_labels)
            names_with_ip = filter(n -> startswith(String(n), "ip"), names(@__MODULE__; all=true))
            max_ip = maximum([getfield(@__MODULE__, n) for n in names_with_ip])
            labels = Vector{String}(undef, max_ip)
            for n in names_with_ip
                  idx = getfield(@__MODULE__, n)
                  labels[idx] = replace(String(n), "ip" => "")
            end
            const ip_labels = labels
      end
      
      const iptotal = 17 
      const ipclimbn   = ipclimb5      # last climb   point
      const ipcruisen  = ipcruise2     # last cruise  point
      const ipdescentn = ipdescent5    # last descent point

#---- keywords for ij sequence parameters 
#=
      character*2 cplab(iptotal)
#---- printout labels for mission points
      data cplab /
      'ST'
      'RO'
      'TO'
      'CB'
      'B1'
      'B2'
      'B3'
      'B4'
      'B5'
      'C1'
      'C2'
      'D1'
      'D2'
      'D3'
      'D4'
      'D5'
      'TE'  /
=#
;