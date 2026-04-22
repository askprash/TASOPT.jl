"""
Enums for turbofan engine calculation options and station identifiers.
"""

using EnumX

# ---------------------------------------------------------------------------
# EngineStation — compile-time identifiers for every named turbofan station
# ---------------------------------------------------------------------------

"""
    EngineStation

Enumeration of every named flow station in the TASOPT turbofan.  Each member
is a compile-time constant that can be used as a typed sentinel to index
engine-state containers, eliminating bare integer indices.

SAE ARP755 station numbering convention:

| Member            | ARP755 # | Description                                  |
|:------------------|:--------:|:---------------------------------------------|
| `Freestream`      | 0        | Far-field ambient                            |
| `FanFaceOuter`    | 12       | Fan face outside casing boundary layers      |
| `FanFaceLPC`      | 2a       | Fan face over LPC (low-pressure compressor)  |
| `PreCoolerOut`    | 2ac      | Pre-cooler outlet / LPC inlet                |
| `FanFaceFan`      | 2        | Fan face over fan stream                     |
| `FanExit`         | 13       | Fan exit / pre-cooler inlet                  |
| `LPCExit`         | 25       | LPC exit / inter-cooler inlet                |
| `InterCoolerOut`  | 25c      | Inter-cooler outlet / HPC inlet              |
| `HPCExit`         | 3        | HPC (high-pressure compressor) exit          |
| `CombustorExit`   | 4        | Combustor exit before cooling air addition   |
| `CoolMixInlet`    | 4a       | Start-of-mixing / cooling-flow outlet        |
| `TurbineInlet`    | 41       | Turbine inlet after cooling air addition     |
| `HPTExit`         | 45       | HPT (high-pressure turbine) exit / LPT inlet |
| `LPTExit`         | 5        | LPT (low-pressure turbine) exit              |
| `RegenCoolerOut`  | 5c       | Regenerative cooler outlet                   |
| `CoreNozzle`      | 8        | Core nozzle throat                           |
| `CoreNozzleExit`  | 9        | Core flow downstream of nozzle               |
| `FanNozzle`       | 18       | Fan nozzle throat                            |
| `FanNozzleExit`   | 19       | Fan duct exit downstream of fan nozzle       |
| `OfftakeDisch`    | 25off    | Offtake air discharge point                  |
"""
@enumx EngineStation begin
    Freestream       # station 0:     far-field ambient
    FanFaceOuter     # station 12:    fan face outside casing boundary layers
    FanFaceLPC       # station 2a:    fan face over LPC portion
    PreCoolerOut     # station 2ac:   pre-cooler outlet, LPC inlet
    FanFaceFan       # station 2:     fan face over fan stream
    FanExit          # station 13:    fan exit (off-design: precooler inlet)
    LPCExit          # station 25:    LPC exit (off-design: intercooler inlet)
    InterCoolerOut   # station 25c:   inter-cooler outlet, HPC inlet
    HPCExit          # station 3:     HPC discharge, combustor inlet
    CombustorExit    # station 4:     combustor exit before cooling air addition
    CoolMixInlet     # station 4a:    start-of-mixing / cooling-flow outlet
    TurbineInlet     # station 41:    turbine inlet after cooling air addition
    HPTExit          # station 45:    HPT exit, LPT inlet
    LPTExit          # station 5:     LPT exit (off-design: regen cooler inlet)
    RegenCoolerOut   # station 5c:    regenerative cooler outlet
    CoreNozzle       # station 8:     core nozzle throat
    CoreNozzleExit   # station 9:     core flow downstream of nozzle
    FanNozzle        # station 18:    fan nozzle throat
    FanNozzleExit    # station 19:    fan duct exit downstream of fan nozzle
    OfftakeDisch     # station 25off: offtake air discharge point
end

"""
    station_number(s::EngineStation.T) -> String

Return the SAE ARP755 station number string for station `s`.  Non-integer
stations (e.g. `"2ac"`, `"4a"`, `"25off"`) are returned as strings.

# Examples
```julia
station_number(EngineStation.Freestream)   # "0"
station_number(EngineStation.HPCExit)      # "3"
station_number(EngineStation.PreCoolerOut) # "2ac"
station_number(EngineStation.LPCExit)      # "25"
station_number(EngineStation.TurbineInlet) # "41"
```
"""
function station_number(s::EngineStation.T)::String
    s == EngineStation.Freestream      && return "0"
    s == EngineStation.FanFaceOuter    && return "12"
    s == EngineStation.FanFaceLPC      && return "2a"
    s == EngineStation.PreCoolerOut    && return "2ac"
    s == EngineStation.FanFaceFan      && return "2"
    s == EngineStation.FanExit         && return "13"
    s == EngineStation.LPCExit         && return "25"
    s == EngineStation.InterCoolerOut  && return "25c"
    s == EngineStation.HPCExit         && return "3"
    s == EngineStation.CombustorExit   && return "4"
    s == EngineStation.CoolMixInlet    && return "4a"
    s == EngineStation.TurbineInlet    && return "41"
    s == EngineStation.HPTExit         && return "45"
    s == EngineStation.LPTExit         && return "5"
    s == EngineStation.RegenCoolerOut  && return "5c"
    s == EngineStation.CoreNozzle      && return "8"
    s == EngineStation.CoreNozzleExit  && return "9"
    s == EngineStation.FanNozzle       && return "18"
    s == EngineStation.FanNozzleExit   && return "19"
    s == EngineStation.OfftakeDisch    && return "25off"
    error("Unknown EngineStation variant: $s")
end

"""
    station_description(s::EngineStation.T) -> String

Return a short human-readable description of station `s`.
"""
function station_description(s::EngineStation.T)::String
    s == EngineStation.Freestream      && return "freestream (ambient)"
    s == EngineStation.FanFaceOuter    && return "fan face outside casing boundary layers"
    s == EngineStation.FanFaceLPC      && return "fan face over LPC portion"
    s == EngineStation.PreCoolerOut    && return "pre-cooler outlet / LPC inlet"
    s == EngineStation.FanFaceFan      && return "fan face over fan stream"
    s == EngineStation.FanExit         && return "fan exit / pre-cooler inlet"
    s == EngineStation.LPCExit         && return "LPC exit / inter-cooler inlet"
    s == EngineStation.InterCoolerOut  && return "inter-cooler outlet / HPC inlet"
    s == EngineStation.HPCExit         && return "HPC discharge / combustor inlet"
    s == EngineStation.CombustorExit   && return "combustor exit before cooling"
    s == EngineStation.CoolMixInlet    && return "start-of-mixing / cooling-flow outlet"
    s == EngineStation.TurbineInlet    && return "turbine inlet after cooling air"
    s == EngineStation.HPTExit         && return "HPT exit / LPT inlet"
    s == EngineStation.LPTExit         && return "LPT exit / regen cooler inlet"
    s == EngineStation.RegenCoolerOut  && return "regenerative cooler outlet"
    s == EngineStation.CoreNozzle      && return "core nozzle throat"
    s == EngineStation.CoreNozzleExit  && return "core nozzle exit / downstream"
    s == EngineStation.FanNozzle       && return "fan nozzle throat"
    s == EngineStation.FanNozzleExit   && return "fan nozzle exit / downstream"
    s == EngineStation.OfftakeDisch    && return "offtake air discharge"
    error("Unknown EngineStation variant: $s")
end

"""
    CalcMode

Selects which turbofan calculation to perform.

- `Sizing`: on-design sizing (`tfsize!`)
- `FixedTt4OffDes`: off-design at fixed turbine inlet temperature
- `FixedFeOffDes`: off-design at fixed net thrust
"""
@enumx CalcMode Sizing FixedTt4OffDes FixedFeOffDes

"""
    CoolingOpt

Selects the turbine cooling model.

- `NoCooling`: no cooling air, station 4.1 == station 4
- `FixedCoolingFlowRatio`: cooling bypass ratios `epsrow` are inputs; `Tmrow` is computed
- `FixedTmetal`: metal temperatures `Tmrow` are inputs; `epsrow` is computed
"""
@enumx CoolingOpt NoCooling FixedCoolingFlowRatio FixedTmetal
