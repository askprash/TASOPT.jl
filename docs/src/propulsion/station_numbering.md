# [Engine Station Numbering](@id station_numbering)

TASOPT uses the SAE ARP755 station numbering convention as the base,
extended with a small number of additional stations for its heat-exchanger
architecture (pre-cooler, inter-cooler, regenerative cooler) and bleed-air
offtake.  This page documents every station defined in the
[`EngineStation`](@ref) enum.

## Station Table

| `EngineStation` member | Station # | ARP755? | Description                                   |
|:-----------------------|:---------:|:-------:|:----------------------------------------------|
| `Freestream`           | 0         | yes     | Far-field ambient                             |
| `FanFaceOuter`         | 12        | yes     | Fan face outside casing boundary layers       |
| `FanFaceLPC`           | 2a        | yes     | Fan face over LPC (low-pressure compressor)   |
| `PreCoolerOut`         | 2ac       | **no**  | Pre-cooler outlet / LPC inlet                 |
| `FanFaceFan`           | 2         | yes     | Fan face over fan stream                      |
| `FanExit`              | 13        | yes     | Fan exit / pre-cooler inlet                   |
| `LPCExit`              | 25        | yes     | LPC exit / inter-cooler inlet                 |
| `InterCoolerOut`       | 25c       | **no**  | Inter-cooler outlet / HPC inlet               |
| `HPCExit`              | 3         | yes     | HPC (high-pressure compressor) exit           |
| `CombustorExit`        | 4         | yes     | Combustor exit before cooling air addition    |
| `CoolMixInlet`         | 4a        | yes     | Start-of-mixing / cooling-flow outlet         |
| `TurbineInlet`         | 41        | yes     | Turbine inlet after cooling air addition      |
| `HPTExit`              | 45        | yes     | HPT (high-pressure turbine) exit / LPT inlet  |
| `LPTExit`              | 5         | yes     | LPT (low-pressure turbine) exit               |
| `RegenCoolerOut`       | 5c        | **no**  | Regenerative cooler outlet                    |
| `CoreNozzle`           | 8         | yes     | Core nozzle throat                            |
| `CoreNozzleExit`       | 9         | yes     | Core flow downstream of nozzle                |
| `FanNozzle`            | 18        | yes     | Fan nozzle throat                             |
| `FanNozzleExit`        | 19        | yes     | Fan duct exit downstream of fan nozzle        |
| `OfftakeDisch`         | 25off     | **no**  | Offtake air discharge point                   |

## TASOPT Extensions

Four stations are **not** defined in ARP755 and are unique to TASOPT:

| Station # | Member          | Purpose                                                     |
|:---------:|:----------------|:------------------------------------------------------------|
| `2ac`     | `PreCoolerOut`  | Exit of the pre-cooler heat exchanger, between fan and LPC  |
| `25c`     | `InterCoolerOut`| Exit of the inter-cooler heat exchanger, between LPC and HPC|
| `5c`      | `RegenCoolerOut`| Exit of the regenerative cooler, downstream of LPT          |
| `25off`   | `OfftakeDisch`  | Bleed-air offtake discharge point                           |

These are numbered by inserting a suffix letter (`c` for cooler, `off` for
offtake) after the nearest upstream ARP755 station number.

## Reference

- SAE ARP755 — *Aircraft Propulsion System Performance Station Designation and
  Nomenclature* (current revision).
- Drela, M., *TASOPT 2.16 documentation*, MIT, 2010 — original station
  definitions for the heat-exchanger architecture.
