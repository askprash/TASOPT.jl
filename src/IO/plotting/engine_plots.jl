"""
engine_plots.jl — Standard engineering plots for engine sweep results.

Provides two functions exported from the TASOPT module:

- [`plot_engine_station_profiles`](@ref): total-temperature and total-pressure
  profiles along the gas-path stations for every mission point.
- [`plot_engine_performance`](@ref): top-level engine performance quantities
  (thrust, TSFC, BPR, core mass flow, and optionally spool speeds) versus
  mission point.

These functions accept either a `Mission{T}` (preferred) or the deprecated
[`engine.SweepResult`](@ref) returned by [`engine.run_engine_sweep`](@ref),
and return a `Plots.Plot` object that can be displayed or saved with `savefig`.
"""

# Gas-path station ordering for station-profile plots.
# Only stations that have a total-temperature entry in `pare` are included;
# HX intermediates (19c, 25c, 4a, 49c) and static-only stations (6, 8) are
# omitted because they are zero for the default turbofan.
const _PROFILE_STATIONS = (
    ("0",  :st0),   # Freestream
    ("18", :st18),  # FanFaceOuter
    ("19", :st19),  # FanFaceLPC
    ("2",  :st2),   # FanFaceFan / LPC inlet
    ("21", :st21),  # FanExit
    ("25", :st25),  # LPCExit / HPC inlet
    ("3",  :st3),   # HPCExit / Combustor inlet
    ("4",  :st4),   # CombustorExit / HPT inlet (Tt4 target)
    ("41", :st41),  # TurbineInlet (post-cooling mix)
    ("45", :st45),  # HPTExit / LPT inlet
    ("49", :st49),  # LPTExit / Nozzle inlet
    ("5",  :st5),   # CoreNozzle throat
    ("7",  :st7),   # FanNozzle throat
)

"""
    plot_engine_station_profiles(mission::Mission, ip_range; kwargs...) -> Plots.Plot

Return a two-panel figure with total-temperature (K, top panel) and
total-pressure (Pa, bottom panel, log-scale) profiles along the engine
gas path for every mission point in `ip_range`.

One line is drawn per mission point; colours follow the default Plots.jl
palette.  Only stations with a meaningful total-state entry in `pare` are
shown (freestream through fan nozzle, 13 stations total).

`kwargs` are forwarded to the outer `Plots.plot` call (e.g. `size`, `dpi`).

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)
mission = TASOPT.engine.run_engine_sweep(ac)
p = TASOPT.plot_engine_station_profiles(mission, ipstatic:ipdescentn)
savefig(p, "station_profiles.png")
```
"""
function plot_engine_station_profiles(mission::Mission, ip_range; kwargs...)
    ns     = length(_PROFILE_STATIONS)
    xlbls  = [t[1] for t in _PROFILE_STATIONS]
    flds   = [t[2] for t in _PROFILE_STATIONS]
    xrange = 1:ns

    p_Tt = Plots.plot(;
        xlabel    = "Station",
        ylabel    = "Tₜ  (K)",
        title     = "Total-temperature profile",
        xticks    = (xrange, xlbls),
        legend    = :topright,
        framestyle = :box,
    )
    p_pt = Plots.plot(;
        xlabel    = "Station",
        ylabel    = "pₜ  (Pa)",
        title     = "Total-pressure profile",
        xticks    = (xrange, xlbls),
        yscale    = :log10,
        legend    = :topright,
        framestyle = :box,
    )

    for ip in ip_range
        eng   = mission.points[ip].engine
        lbl   = (ip <= length(ip_labels)) ? ip_labels[ip] : string(ip)
        Tt_v  = [getfield(eng, f).Tt for f in flds]
        pt_v  = [getfield(eng, f).pt for f in flds]
        Plots.plot!(p_Tt, xrange, Tt_v; label=lbl, markershape=:circle, ms=3)
        Plots.plot!(p_pt, xrange, pt_v; label=lbl, markershape=:circle, ms=3)
    end

    Plots.plot(p_Tt, p_pt; layout=(2, 1), kwargs...)
end

"""
    plot_engine_station_profiles(result::SweepResult; kwargs...) -> Plots.Plot

!!! warning "Deprecated"
    This overload is deprecated.  Use
    `plot_engine_station_profiles(mission, ip_range)` with the `Mission{T}`
    returned by [`engine.run_engine_sweep`](@ref).
"""
function plot_engine_station_profiles(result::engine.SweepResult; kwargs...)
    ns     = length(_PROFILE_STATIONS)
    xlbls  = [t[1] for t in _PROFILE_STATIONS]
    flds   = [t[2] for t in _PROFILE_STATIONS]
    xrange = 1:ns
    n      = length(result.ip_indices)

    p_Tt = Plots.plot(;
        xlabel    = "Station",
        ylabel    = "Tₜ  (K)",
        title     = "Total-temperature profile",
        xticks    = (xrange, xlbls),
        legend    = :topright,
        framestyle = :box,
    )
    p_pt = Plots.plot(;
        xlabel    = "Station",
        ylabel    = "pₜ  (Pa)",
        title     = "Total-pressure profile",
        xticks    = (xrange, xlbls),
        yscale    = :log10,
        legend    = :topright,
        framestyle = :box,
    )

    for k in 1:n
        eng   = result.engines[k]
        lbl   = result.ip_labels[k]
        Tt_v  = [getfield(eng, f).Tt for f in flds]
        pt_v  = [getfield(eng, f).pt for f in flds]
        Plots.plot!(p_Tt, xrange, Tt_v; label=lbl, markershape=:circle, ms=3)
        Plots.plot!(p_pt, xrange, pt_v; label=lbl, markershape=:circle, ms=3)
    end

    Plots.plot(p_Tt, p_pt; layout=(2, 1), kwargs...)
end

"""
    plot_engine_performance(mission::Mission, ip_range; ac=nothing, imission=1,
                            kwargs...) -> Plots.Plot

Return a multi-panel figure of top-level engine performance quantities versus
mission point.

**Always included** (four panels):

| Panel | Quantity | Unit |
|:------|:---------|:-----|
| 1 | Net thrust Fₑ | N |
| 2 | TSFC | mg/N/s |
| 3 | Bypass ratio BPR | — |
| 4 | Core mass flow ṁ_core | kg/s |

**Optional spool-speed panels** (three additional panels, requires `ac`):

| Panel | Quantity | Unit |
|:------|:---------|:-----|
| 5 | Fan corrected speed N̄_f | — |
| 6 | LPC corrected speed N̄_lc | — |
| 7 | HPC corrected speed N̄_hc | — |

`imission` selects the mission index used for spool-speed lookup (default 1).

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)
mission = TASOPT.engine.run_engine_sweep(ac)

# Basic four-panel plot (no spool speeds)
p_basic = TASOPT.plot_engine_performance(mission, ipstatic:ipdescentn)

# Full seven-panel plot (with spool speeds)
p_full = TASOPT.plot_engine_performance(mission, ipstatic:ipdescentn; ac=ac)
savefig(p_full, "engine_performance.png")
```
"""
function plot_engine_performance(mission::Mission, ip_range;
                                 ac=nothing, imission::Int=1, kwargs...)
    ips    = collect(Int, ip_range)
    n      = length(ips)
    xrange = 1:n

    lbls = [(ip <= length(ip_labels)) ? ip_labels[ip] : string(ip) for ip in ips]
    xtk  = (xrange, lbls)

    Fe_v    = [mission.points[ip].engine.Fe      for ip in ips]
    TSFC_v  = [mission.points[ip].engine.TSFC    for ip in ips]
    BPR_v   = [mission.points[ip].engine.BPR     for ip in ips]
    mcore_v = [mission.points[ip].engine.st2.mdot for ip in ips]
    tsfc_mg = TSFC_v .* 1f6

    p_Fe   = Plots.bar(xrange, Fe_v;
                       ylabel="Fe (N)", title="Net thrust",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_TSFC = Plots.bar(xrange, tsfc_mg;
                       ylabel="TSFC (mg/N/s)", title="Fuel consumption",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_BPR  = Plots.bar(xrange, BPR_v;
                       ylabel="BPR", title="Bypass ratio",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_mc   = Plots.bar(xrange, mcore_v;
                       ylabel="mcore (kg/s)", title="Core mass flow",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)

    if ac === nothing
        return Plots.plot(p_Fe, p_TSFC, p_BPR, p_mc;
                          layout=(2, 2), kwargs...)
    end

    Nbf_v  = [ac.pare[ieNf,   ip, imission] for ip in ips]
    Nblc_v = [ac.pare[ieNblc, ip, imission] for ip in ips]
    Nbhc_v = [ac.pare[ieNbhc, ip, imission] for ip in ips]

    p_Nbf  = Plots.bar(xrange, Nbf_v;
                       ylabel="Nbf", title="Fan spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_Nblc = Plots.bar(xrange, Nblc_v;
                       ylabel="Nblc", title="LPC spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_Nbhc = Plots.bar(xrange, Nbhc_v;
                       ylabel="Nbhc", title="HPC spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)

    Plots.plot(p_Fe, p_TSFC, p_BPR, p_mc, p_Nbf, p_Nblc, p_Nbhc;
               layout=(4, 2), kwargs...)
end

"""
    plot_engine_performance(result::SweepResult; ac=nothing, imission=1,
                            kwargs...) -> Plots.Plot

!!! warning "Deprecated"
    This overload is deprecated.  Use
    `plot_engine_performance(mission, ip_range)` with the `Mission{T}`
    returned by [`engine.run_engine_sweep`](@ref).
"""
function plot_engine_performance(result::engine.SweepResult;
                                 ac=nothing, imission::Int=1, kwargs...)
    n      = length(result.ip_indices)
    xrange = 1:n

    xtk = (xrange, result.ip_labels)

    tsfc_mg = result.TSFC .* 1f6

    p_Fe   = Plots.bar(xrange, result.Fe;
                       ylabel="Fe (N)", title="Net thrust",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_TSFC = Plots.bar(xrange, tsfc_mg;
                       ylabel="TSFC (mg/N/s)", title="Fuel consumption",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_BPR  = Plots.bar(xrange, result.BPR;
                       ylabel="BPR", title="Bypass ratio",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_mc   = Plots.bar(xrange, result.mcore;
                       ylabel="mcore (kg/s)", title="Core mass flow",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)

    if ac === nothing
        return Plots.plot(p_Fe, p_TSFC, p_BPR, p_mc;
                          layout=(2, 2), kwargs...)
    end

    Nbf_v  = [ac.pare[ieNf,   ip, imission] for ip in result.ip_indices]
    Nblc_v = [ac.pare[ieNblc, ip, imission] for ip in result.ip_indices]
    Nbhc_v = [ac.pare[ieNbhc, ip, imission] for ip in result.ip_indices]

    p_Nbf  = Plots.bar(xrange, Nbf_v;
                       ylabel="Nbf", title="Fan spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_Nblc = Plots.bar(xrange, Nblc_v;
                       ylabel="Nblc", title="LPC spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)
    p_Nbhc = Plots.bar(xrange, Nbhc_v;
                       ylabel="Nbhc", title="HPC spool speed (corr.)",
                       xticks=xtk, xrotation=45, xlabel="Mission point", framestyle=:box)

    Plots.plot(p_Fe, p_TSFC, p_BPR, p_mc, p_Nbf, p_Nblc, p_Nbhc;
               layout=(4, 2), kwargs...)
end
