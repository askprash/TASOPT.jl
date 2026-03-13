### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of
# Pluto, the following 'mock version' of @bind gives bound variables a default value
# (instead of throwing an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a0000001-0000-0000-0000-000000000001
begin
    import Pkg
    # Activate the TASOPT project (notebook lives in example/ subfolder)
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.add(["PlutoUI", "StatsPlots"])
    Pkg.instantiate()
end

# ╔═╡ a0000002-0000-0000-0000-000000000002
begin
    using TASOPT
    using Plots
    using StatsPlots
    using PlutoUI
    using Printf
    include(TASOPT.__TASOPTindices__)
    TableOfContents(title="✈ Comparison Dashboard")
end

# ╔═╡ a0000003-0000-0000-0000-000000000003
md"""
# ✈ TASOPT Aircraft Comparison Dashboard

Compare two sized TASOPT aircraft designs side-by-side. Load saved `.jld2` files
(created with `TASOPT.quicksave_aircraft(ac, "myfile.jld2")`) using the inputs below.

---
"""

# ╔═╡ a0000004-0000-0000-0000-000000000004
md"""
## 📂 Aircraft File Selection

| | Path |
|:--|:--|
| **Aircraft 1** | $(@bind ac1_path TextField(80; default=joinpath(TASOPT.__TASOPTroot__, "IO/default_quicksave_aircraft.jld2"))) |
| **Aircraft 2** | $(@bind ac2_path TextField(80; default=joinpath(TASOPT.__TASOPTroot__, "IO/default_quicksave_aircraft.jld2"))) |

> **Tip:** Save a sized aircraft with `TASOPT.quicksave_aircraft(ac, "path/to/file.jld2")`, then paste the path above.
"""

# ╔═╡ a0000005-0000-0000-0000-000000000005
begin
    function try_load(path)
        isfile(path) || return (nothing, "File not found: $path")
        try
            ac = TASOPT.quickload_aircraft(path)
            ac.is_sized[1] || return (nothing, "Aircraft at $path is not sized — run size_aircraft! first")
            return (ac, nothing)
        catch e
            return (nothing, "Error loading $path: $e")
        end
    end

    ac1, err1 = try_load(ac1_path)
    ac2, err2 = try_load(ac2_path)

    if !isnothing(err1)
        @warn err1
    end
    if !isnothing(err2)
        @warn err2
    end
end;

# ╔═╡ a0000006-0000-0000-0000-000000000006
# Status banner
if isnothing(ac1) || isnothing(ac2)
    md"""
    !!! danger "Aircraft not loaded"
        $(isnothing(ac1) ? "**Aircraft 1:** $(err1)" : "")
        $(isnothing(ac2) ? "**Aircraft 2:** $(err2)" : "")
    """
else
    md"""
    !!! success "Both aircraft loaded successfully"
        **Aircraft 1:** $(ac1.name)
        **Aircraft 2:** $(ac2.name)
    """
end

# ╔═╡ a0000007-0000-0000-0000-000000000007
md"""
---
## 📊 Overview Comparison
"""

# ╔═╡ a0000008-0000-0000-0000-000000000008
"""
    fmt_delta(v1, v2; pct=true)

Return a formatted Δ string showing the difference of v2 relative to v1.
"""
function fmt_delta(v1, v2; pct=true)
    v1 == 0.0 && return "—"
    d = (v2 - v1) / abs(v1) * 100
    color = d > 0 ? "green" : (d < 0 ? "red" : "gray")
    sign  = d > 0 ? "+" : ""
    return """<span style="color:$(color)">$(sign)$(round(d, digits=1))%</span>"""
end

# ╔═╡ a0000009-0000-0000-0000-000000000009
# Overview comparison table
if !isnothing(ac1) && !isnothing(ac2)
    # Collect metrics
    lbf_to_N = 4.44822

    # Mass/weight metrics (kN → t via /g)
    WMTO1  = ac1.parg[igWMTO];   WMTO2  = ac2.parg[igWMTO]
    Wfuel1 = ac1.parg[igWfuel];  Wfuel2 = ac2.parg[igWfuel]
    Wpay1  = ac1.parg[igWpay];   Wpay2  = ac2.parg[igWpay]
    Wempty1 = WMTO1 - Wfuel1 - Wpay1
    Wempty2 = WMTO2 - Wfuel2 - Wpay2

    # Mission performance
    PFEI1  = ac1.parm[imPFEI, 1];   PFEI2  = ac2.parm[imPFEI, 1]
    Range1 = ac1.parm[imRange, 1];  Range2 = ac2.parm[imRange, 1]

    # Aerodynamics at cruise
    CL1 = ac1.para[iaCL, ipcruise1, 1];  CL2 = ac2.para[iaCL, ipcruise1, 1]
    CD1 = ac1.para[iaCD, ipcruise1, 1];  CD2 = ac2.para[iaCD, ipcruise1, 1]
    LoD1 = CL1/CD1;  LoD2 = CL2/CD2
    M1 = ac1.para[iaMach, ipcruise1, 1]; M2 = ac2.para[iaMach, ipcruise1, 1]

    # Geometry
    b1 = ac1.wing.layout.span;       b2 = ac2.wing.layout.span
    AR1 = ac1.wing.layout.AR;        AR2 = ac2.wing.layout.AR
    S1 = ac1.wing.layout.S;          S2 = ac2.wing.layout.S
    sweep1 = ac1.wing.layout.sweep;  sweep2 = ac2.wing.layout.sweep
    Rfuse1 = ac1.fuselage.layout.radius; Rfuse2 = ac2.fuselage.layout.radius
    Lfuse1 = ac1.fuselage.layout.x_end - ac1.fuselage.layout.x_nose
    Lfuse2 = ac2.fuselage.layout.x_end - ac2.fuselage.layout.x_nose

    # Engine
    dfan1 = ac1.parg[igdfan];  dfan2 = ac2.parg[igdfan]
    neng1 = ac1.parg[igneng];  neng2 = ac2.parg[igneng]
    BPR1 = ac1.pare[ieBPR, ipcruise1, 1]; BPR2 = ac2.pare[ieBPR, ipcruise1, 1]
    OPR1 = ac1.pare[ieOPR, ipcruise1, 1]; OPR2 = ac2.pare[ieOPR, ipcruise1, 1]
    TSFC1 = ac1.pare[ieTSFC, ipcruise1, 1]; TSFC2 = ac2.pare[ieTSFC, ipcruise1, 1]

    rows = [
        ("**METRIC**",                         "**$(ac1.name)**",                          "**$(ac2.name)**",                          "**Δ (2 vs 1)**"),
        ("*Performance*",                       "",                                          "",                                          ""),
        ("PFEI [kJ/kg/km]",                    @sprintf("%.3f", PFEI1),                    @sprintf("%.3f", PFEI2),                    fmt_delta(PFEI1, PFEI2)),
        ("Design Range [km]",                  @sprintf("%.0f", Range1/1e3),               @sprintf("%.0f", Range2/1e3),               fmt_delta(Range1, Range2)),
        ("Cruise L/D",                         @sprintf("%.2f", LoD1),                     @sprintf("%.2f", LoD2),                     fmt_delta(LoD1, LoD2)),
        ("Cruise Mach",                        @sprintf("%.3f", M1),                       @sprintf("%.3f", M2),                       fmt_delta(M1, M2)),
        ("*Weights*",                           "",                                          "",                                          ""),
        ("MTOW [t]",                           @sprintf("%.1f", WMTO1/9.81/1e3),           @sprintf("%.1f", WMTO2/9.81/1e3),           fmt_delta(WMTO1, WMTO2)),
        ("OEW [t]",                            @sprintf("%.1f", Wempty1/9.81/1e3),         @sprintf("%.1f", Wempty2/9.81/1e3),         fmt_delta(Wempty1, Wempty2)),
        ("Fuel [t]",                           @sprintf("%.1f", Wfuel1/9.81/1e3),          @sprintf("%.1f", Wfuel2/9.81/1e3),          fmt_delta(Wfuel1, Wfuel2)),
        ("Payload [t]",                        @sprintf("%.1f", Wpay1/9.81/1e3),           @sprintf("%.1f", Wpay2/9.81/1e3),           fmt_delta(Wpay1, Wpay2)),
        ("*Wing Geometry*",                     "",                                          "",                                          ""),
        ("Span [m]",                           @sprintf("%.1f", b1),                       @sprintf("%.1f", b2),                       fmt_delta(b1, b2)),
        ("Aspect Ratio",                       @sprintf("%.2f", AR1),                      @sprintf("%.2f", AR2),                      fmt_delta(AR1, AR2)),
        ("Wing Area [m²]",                     @sprintf("%.1f", S1),                       @sprintf("%.1f", S2),                       fmt_delta(S1, S2)),
        ("Sweep [°]",                          @sprintf("%.1f", sweep1),                   @sprintf("%.1f", sweep2),                   fmt_delta(sweep1, sweep2)),
        ("*Fuselage*",                          "",                                          "",                                          ""),
        ("Fuselage Length [m]",                @sprintf("%.1f", Lfuse1),                   @sprintf("%.1f", Lfuse2),                   fmt_delta(Lfuse1, Lfuse2)),
        ("Fuselage Radius [m]",                @sprintf("%.2f", Rfuse1),                   @sprintf("%.2f", Rfuse2),                   fmt_delta(Rfuse1, Rfuse2)),
        ("*Engine*",                            "",                                          "",                                          ""),
        ("No. Engines",                        @sprintf("%.0f", neng1),                    @sprintf("%.0f", neng2),                    fmt_delta(neng1, neng2)),
        ("Fan Diameter [m]",                   @sprintf("%.2f", dfan1),                    @sprintf("%.2f", dfan2),                    fmt_delta(dfan1, dfan2)),
        ("BPR",                                @sprintf("%.1f", BPR1),                     @sprintf("%.1f", BPR2),                     fmt_delta(BPR1, BPR2)),
        ("OPR",                                @sprintf("%.1f", OPR1),                     @sprintf("%.1f", OPR2),                     fmt_delta(OPR1, OPR2)),
        ("TSFC [mg/N/s]",                      @sprintf("%.4f", TSFC1*1e6),               @sprintf("%.4f", TSFC2*1e6),               fmt_delta(TSFC1, TSFC2)),
    ]

    # Build HTML table with striped rows
    html_str = """
    <style>
      .actable { border-collapse: collapse; width: 100%; font-family: monospace; font-size: 0.9em; }
      .actable th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }
      .actable td { padding: 6px 12px; border-bottom: 1px solid #ddd; }
      .actable tr:nth-child(even) { background: #f8f9fa; }
      .actable tr.section-header td { background: #ecf0f1; font-weight: bold; color: #2c3e50; }
      .actable tr:hover { background: #e8f4f8; }
    </style>
    <table class="actable">
    """

    for (i, row) in enumerate(rows)
        if i == 1
            html_str *= "<tr>" * join(["<th>$c</th>" for c in row]) * "</tr>\n"
        elseif isempty(strip(row[2]))  # section header row
            html_str *= """<tr class="section-header">""" * join(["<td>$(replace(c, r"\*" => ""))</td>" for c in row]) * "</tr>\n"
        else
            html_str *= "<tr>" * join(["<td>$c</td>" for c in row]) * "</tr>\n"
        end
    end
    html_str *= "</table>"

    HTML(html_str)
end

# ╔═╡ a0000010-0000-0000-0000-000000000010
md"""
---
## ✏️ Aircraft Geometry (Top View)

Side-by-side stick figures showing the planform geometry of each aircraft.
"""

# ╔═╡ a0000011-0000-0000-0000-000000000011
if !isnothing(ac1) && !isnothing(ac2)
    p1 = TASOPT.stickfig(ac1;
        annotate_text=true, annotate_length=true, annotate_group=false,
        label_fs=10)
    title!(p1, ac1.name)

    p2 = TASOPT.stickfig(ac2;
        annotate_text=true, annotate_length=true, annotate_group=false,
        label_fs=10)
    title!(p2, ac2.name)

    plot(p1, p2;
        layout=(2, 1),
        size=(900, 900),
        dpi=150)
end

# ╔═╡ a0000012-0000-0000-0000-000000000012
md"""
---
## ⚖️ Weight Breakdown
"""

# ╔═╡ a0000013-0000-0000-0000-000000000013
"""
    extract_weights(ac) -> (labels, values_N)

Extract component weights from an aircraft struct [N].
"""
function extract_weights(ac)
    parg = ac.parg
    fuse = ac.fuselage
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    lg = ac.landing_gear

    Wbox  = wing.inboard.webs.weight.W + wing.inboard.caps.weight.W
    Wwing = Wbox * (1 + wing.weight_frac_flap + wing.weight_frac_slat +
                    wing.weight_frac_ailerons + wing.weight_frac_leading_trailing_edge +
                    wing.weight_frac_ribs + wing.weight_frac_spoilers +
                    wing.weight_frac_attachments)

    Whpe  = parg[igWMTO] * fuse.HPE_sys.W
    Wadd  = Whpe + lg.nose_gear.weight.W + lg.main_gear.weight.W

    labels = ["Fuselage", "Wing", "H-Tail", "V-Tail", "Engine Sys", "Fuel Tank", "Added"]
    vals   = [fuse.weight, Wwing, htail.weight, vtail.weight,
              parg[igWtesys], parg[igWftank], Wadd]
    return labels, vals
end

# ╔═╡ a0000014-0000-0000-0000-000000000014
if !isnothing(ac1) && !isnothing(ac2)
    wlabels, wvals1 = extract_weights(ac1)
    _,       wvals2 = extract_weights(ac2)

    # Convert to tonnes (÷ 9.81 × 1e-3)
    w1_t = wvals1 ./ 9.81e3
    w2_t = wvals2 ./ 9.81e3

    p_w = groupedbar(
        wlabels,
        [w1_t w2_t];
        bar_position = :dodge,
        bar_width    = 0.7,
        label        = [ac1.name ac2.name],
        color        = [:steelblue :tomato],
        alpha        = 0.85,
        ylabel       = "Weight [t]",
        title        = "Structural Weight Breakdown",
        legend       = :topright,
        size         = (800, 450),
        xtickfontsize = 9,
        dpi          = 150,
        grid         = :y,
        gridcolor    = :lightgray
    )

    # ---- MTOW / OEW / fuel pie charts side by side ----
    WMTO1 = ac1.parg[igWMTO]; WMTO2 = ac2.parg[igWMTO]
    Wf1   = ac1.parg[igWfuel]; Wf2   = ac2.parg[igWfuel]
    Wp1   = ac1.parg[igWpay];  Wp2   = ac2.parg[igWpay]
    We1   = WMTO1 - Wf1 - Wp1;  We2   = WMTO2 - Wf2 - Wp2

    pie_labels = ["OEW", "Fuel", "Payload"]
    pp1 = pie(pie_labels, [We1, Wf1, Wp1] ./ WMTO1;
              title="$(ac1.name)\nMTOW = $(round(WMTO1/9.81e3, digits=1)) t",
              color=[:steelblue :skyblue :lightblue],
              legend=false, dpi=150, titlefontsize=9)

    pp2 = pie(pie_labels, [We2, Wf2, Wp2] ./ WMTO2;
              title="$(ac2.name)\nMTOW = $(round(WMTO2/9.81e3, digits=1)) t",
              color=[:tomato :salmon :mistyrose],
              legend=false, dpi=150, titlefontsize=9)

    plot(p_w, pp1, pp2;
        layout=@layout([a{0.55w} b c]),
        size=(1100, 420))
end

# ╔═╡ a0000015-0000-0000-0000-000000000015
md"""
---
## 🌊 Aerodynamics — Drag Breakdown at Cruise
"""

# ╔═╡ a0000016-0000-0000-0000-000000000016
if !isnothing(ac1) && !isnothing(ac2)
    p_drag1 = TASOPT.plot_drag_breakdown(ac1;
        ip=[ipcruise1], show_fractions=false, show_values=true)
    title!(p_drag1, ac1.name)

    p_drag2 = TASOPT.plot_drag_breakdown(ac2;
        ip=[ipcruise1], show_fractions=false, show_values=true)
    title!(p_drag2, ac2.name)

    # CD component comparison at cruise
    drag_keys  = ["CDi", "CDfuse", "CDwing", "CDhtail", "CDvtail", "CDnace"]
    drag_inds  = [iaCDi, iaCDfuse, iaCDwing, iaCDhtail, iaCDvtail, iaCDnace]

    cd1 = [ac1.para[idx, ipcruise1, 1] for idx in drag_inds]
    cd2 = [ac2.para[idx, ipcruise1, 1] for idx in drag_inds]

    # Multiply by Sref to get drag areas [m²]
    S1 = ac1.wing.layout.S; S2 = ac2.wing.layout.S
    da1 = cd1 .* S1 .* 1e4   # cm²
    da2 = cd2 .* S2 .* 1e4   # cm²

    p_cd_comp = groupedbar(
        drag_keys,
        [da1 da2];
        bar_position = :dodge,
        bar_width    = 0.7,
        label        = [ac1.name ac2.name],
        color        = [:steelblue :tomato],
        alpha        = 0.85,
        ylabel       = "Drag Area CD×S [cm²]",
        title        = "Cruise Drag Component Comparison",
        legend       = :topright,
        size         = (800, 400),
        dpi          = 150,
        grid         = :y,
        gridcolor    = :lightgray
    )

    plot(p_cd_comp, p_drag1, p_drag2;
        layout=@layout([a{0.45h}; b c]),
        size=(900, 750))
end

# ╔═╡ a0000017-0000-0000-0000-000000000017
md"""
---
## 📈 Mission Performance — Flight Profile
"""

# ╔═╡ a0000018-0000-0000-0000-000000000018
"""
    mission_profile_plot(ac; imission=1, label="", c=:steelblue)

Plot altitude, Mach, CL, and CD over mission flight points for one aircraft.
"""
function mission_profile_plot(ac; imission=1, label="", c=:steelblue)
    npts = size(ac.para, 2)
    pts  = 1:npts

    alt  = ac.para[iaalt,  :, imission] ./ 1e3      # km
    mach = ac.para[iaMach, :, imission]
    CL   = ac.para[iaCL,   :, imission]
    CD   = ac.para[iaCD,   :, imission]
    LoD  = CL ./ CD

    p_alt  = plot(pts, alt;  ylabel="Altitude [km]",  label=label, color=c, lw=2, title="Altitude",  legend=:topright)
    p_mach = plot(pts, mach; ylabel="Mach",           label=label, color=c, lw=2, title="Mach",      legend=false)
    p_LoD  = plot(pts, LoD;  ylabel="L/D",            label=label, color=c, lw=2, title="L/D",       legend=false)
    p_CD   = plot(pts, CD .* 1e4; ylabel="CD [×10⁻⁴]", label=label, color=c, lw=2, title="Drag",    legend=false)

    return p_alt, p_mach, p_LoD, p_CD
end

# ╔═╡ a0000019-0000-0000-0000-000000000019
if !isnothing(ac1) && !isnothing(ac2)
    a1, b1_p, c1, d1 = mission_profile_plot(ac1; label=ac1.name, c=:steelblue)
    a2, b2_p, c2, d2 = mission_profile_plot(ac2; label=ac2.name, c=:tomato)

    # Overlay both aircraft on each panel
    plot!(a1, ac2.para[iaalt,  :, 1] ./ 1e3; label=ac2.name, color=:tomato, lw=2)
    plot!(b1_p, ac2.para[iaMach, :, 1]; label="", color=:tomato, lw=2)
    plot!(c1, ac2.para[iaCL, :, 1] ./ ac2.para[iaCD, :, 1]; label="", color=:tomato, lw=2)
    plot!(d1, ac2.para[iaCD, :, 1] .* 1e4; label="", color=:tomato, lw=2)

    xlabel!.(  [a1, b1_p, c1, d1], "Flight Point")

    plot(a1, b1_p, c1, d1;
        layout=(2, 2),
        size=(900, 550),
        dpi=150,
        leftmargin=5Plots.mm,
        bottommargin=4Plots.mm)
end

# ╔═╡ a0000020-0000-0000-0000-000000000020
md"""
---
## 🔧 Engine Parameters at Cruise
"""

# ╔═╡ a0000021-0000-0000-0000-000000000021
if !isnothing(ac1) && !isnothing(ac2)
    eng_labels = ["BPR", "OPR", "TSFC\n[mg/N/s]", "Fan Dia\n[m]", "Neng"]
    eng_inds   = [ieBPR, ieOPR, ieTSFC, 0, 0]   # 0 = from parg

    function get_eng_vals(ac)
        return [
            ac.pare[ieBPR,  ipcruise1, 1],
            ac.pare[ieOPR,  ipcruise1, 1],
            ac.pare[ieTSFC, ipcruise1, 1] * 1e6,  # convert to mg/N/s
            ac.parg[igdfan],
            ac.parg[igneng]
        ]
    end

    ev1 = get_eng_vals(ac1)
    ev2 = get_eng_vals(ac2)

    # Normalise each parameter so they can be shown on one bar chart
    ev_norm1 = ev1 ./ max.(ev1, ev2)
    ev_norm2 = ev2 ./ max.(ev1, ev2)

    p_eng_norm = groupedbar(
        eng_labels,
        [ev_norm1 ev_norm2];
        bar_position=:dodge, bar_width=0.7,
        label=[ac1.name ac2.name],
        color=[:steelblue :tomato], alpha=0.85,
        ylabel="Normalised value",
        title="Engine Parameters (normalised to max)",
        legend=:bottomright,
        ylims=(0, 1.15),
        size=(700, 380),
        dpi=150,
        grid=:y, gridcolor=:lightgray
    )

    # Absolute value bar charts for BPR and OPR
    p_bpr = bar(
        [ac1.name, ac2.name],
        [ev1[1], ev2[1]];
        color=[:steelblue, :tomato], alpha=0.85,
        ylabel="BPR", title="Bypass Ratio",
        legend=false, dpi=150, ylims=(0, max(ev1[1], ev2[1]) * 1.2)
    )
    p_opr = bar(
        [ac1.name, ac2.name],
        [ev1[2], ev2[2]];
        color=[:steelblue, :tomato], alpha=0.85,
        ylabel="OPR", title="Overall Pressure Ratio",
        legend=false, dpi=150, ylims=(0, max(ev1[2], ev2[2]) * 1.2)
    )
    p_tsfc = bar(
        [ac1.name, ac2.name],
        [ev1[3], ev2[3]];
        color=[:steelblue, :tomato], alpha=0.85,
        ylabel="TSFC [mg/N/s]", title="Thrust-Specific Fuel Consumption",
        legend=false, dpi=150, ylims=(0, max(ev1[3], ev2[3]) * 1.2)
    )

    plot(p_eng_norm, p_bpr, p_opr, p_tsfc;
        layout=@layout([a{0.5h}; b c d]),
        size=(900, 700), dpi=150,
        leftmargin=5Plots.mm)
end

# ╔═╡ a0000022-0000-0000-0000-000000000022
md"""
---
## 📐 Wing & Fuselage Geometry Detail
"""

# ╔═╡ a0000023-0000-0000-0000-000000000023
if !isnothing(ac1) && !isnothing(ac2)
    # ---- Wing chord distribution (normalised) ----
    function wing_chord_dist(ac)
        w = ac.wing
        bo = w.layout.root_span
        bs = w.layout.break_span
        b  = w.layout.span
        co = w.layout.root_chord
        cs = co * w.inboard.λ
        ct = co * w.outboard.λ

        η = [0, bo/b, bs/b, 1.0]
        c = [co, co, cs, ct]   # chord at each station
        return η, c
    end

    η1, c1_w = wing_chord_dist(ac1)
    η2, c2_w = wing_chord_dist(ac2)

    p_chord = plot(η1, c1_w;
        label=ac1.name, color=:steelblue, lw=2.5, marker=:circle, ms=5,
        xlabel="η = 2y/b",  ylabel="Chord [m]",
        title="Wing Chord Distribution",
        legend=:topright, dpi=150, size=(500, 320))
    plot!(p_chord, η2, c2_w;
        label=ac2.name, color=:tomato, lw=2.5, marker=:circle, ms=5)

    # ---- Fuselage cross-section comparison ----
    function fuse_circle(ac)
        R = ac.fuselage.layout.radius
        wfb = ac.fuselage.layout.bubble_center_y_offset
        θ = LinRange(0, 2π, 200)
        x = R .* cos.(θ)
        y = R .* sin.(θ)
        return x, y, R, wfb
    end

    x1c, y1c, R1, _ = fuse_circle(ac1)
    x2c, y2c, R2, _ = fuse_circle(ac2)

    p_fuse_x = plot(x1c, y1c;
        label="$(ac1.name) (R=$(round(R1, digits=2)) m)",
        color=:steelblue, lw=2, aspect_ratio=:equal,
        xlabel="y [m]", ylabel="z [m]",
        title="Fuselage Cross-Section",
        legend=:topright, dpi=150, size=(350, 350))
    plot!(p_fuse_x, x2c, y2c;
        label="$(ac2.name) (R=$(round(R2, digits=2)) m)",
        color=:tomato, lw=2, linestyle=:dash)

    # ---- Key geometry spider/radar-ish: bar chart ----
    geom_labels = ["Span\n[m]", "Fuse L\n[m]", "Wing S\n[m²]", "Root c\n[m]", "Tip c\n[m]"]
    Lfuse1 = ac1.fuselage.layout.x_end - ac1.fuselage.layout.x_nose
    Lfuse2 = ac2.fuselage.layout.x_end - ac2.fuselage.layout.x_nose

    gv1 = [ac1.wing.layout.span, Lfuse1, ac1.wing.layout.S,
           ac1.wing.layout.root_chord, ac1.wing.layout.root_chord * ac1.wing.outboard.λ]
    gv2 = [ac2.wing.layout.span, Lfuse2, ac2.wing.layout.S,
           ac2.wing.layout.root_chord, ac2.wing.layout.root_chord * ac2.wing.outboard.λ]

    p_geom = groupedbar(
        geom_labels,
        [gv1 gv2];
        bar_position=:dodge, bar_width=0.7,
        label=[ac1.name ac2.name],
        color=[:steelblue :tomato], alpha=0.85,
        ylabel="Value [m or m²]",
        title="Key Geometry Parameters",
        legend=:topright, dpi=150, size=(600, 370),
        grid=:y, gridcolor=:lightgray
    )

    plot(p_chord, p_fuse_x, p_geom;
        layout=@layout([a b; c{0.45h}]),
        size=(900, 700), dpi=150,
        leftmargin=5Plots.mm, bottommargin=4Plots.mm)
end

# ╔═╡ a0000024-0000-0000-0000-000000000024
md"""
---
## 🧮 Stability & Balance
"""

# ╔═╡ a0000025-0000-0000-0000-000000000025
if !isnothing(ac1) && !isnothing(ac2)
    function stability_bar(ac, c; label="")
        xNP   = ac.parg[igxNP]
        xCGfw = ac.parg[igxCGfwd]
        xCGaf = ac.parg[igxCGaft]
        mac   = ac.wing.mean_aero_chord

        SMfwd = (xNP - xCGfw) / mac * 100  # %MAC
        SMaft = (xNP - xCGaf) / mac * 100

        bar(["SM fwd (%MAC)", "SM aft (%MAC)"],
            [SMfwd, SMaft];
            color=c, alpha=0.8, label=label,
            ylabel="Static Margin [%MAC]",
            title="$(ac.name)\nNP=$(round(xNP, digits=1)) m  CG fwd=$(round(xCGfw,digits=1)) m  aft=$(round(xCGaf,digits=1)) m",
            dpi=150, titlefontsize=8, legend=false,
            ylims=(0, max(SMfwd, SMaft)*1.4))
    end

    ps1 = stability_bar(ac1, :steelblue)
    ps2 = stability_bar(ac2, :tomato)

    plot(ps1, ps2; layout=(1, 2), size=(750, 380), dpi=150,
         leftmargin=5Plots.mm)
end

# ╔═╡ a0000026-0000-0000-0000-000000000026
md"""
---
*Generated by the TASOPT Aircraft Comparison Dashboard — powered by [Pluto.jl](https://plutojl.org)*
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
# PLUTO PACKAGE MANAGER — do not edit this cell
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
PlutoUI = "0.7"
StatsPlots = "0.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = ""  # populated automatically by Pluto

# ╔═╡ Cell Order:
# ╠═a0000001-0000-0000-0000-000000000001
# ╟─a0000002-0000-0000-0000-000000000002
# ╟─a0000003-0000-0000-0000-000000000003
# ╟─a0000004-0000-0000-0000-000000000004
# ╠═a0000005-0000-0000-0000-000000000005
# ╟─a0000006-0000-0000-0000-000000000006
# ╟─a0000007-0000-0000-0000-000000000007
# ╟─a0000009-0000-0000-0000-000000000009
# ╠═a0000008-0000-0000-0000-000000000008
# ╟─a0000010-0000-0000-0000-000000000010
# ╠═a0000011-0000-0000-0000-000000000011
# ╟─a0000012-0000-0000-0000-000000000012
# ╠═a0000013-0000-0000-0000-000000000013
# ╠═a0000014-0000-0000-0000-000000000014
# ╟─a0000015-0000-0000-0000-000000000015
# ╠═a0000016-0000-0000-0000-000000000016
# ╟─a0000017-0000-0000-0000-000000000017
# ╠═a0000018-0000-0000-0000-000000000018
# ╠═a0000019-0000-0000-0000-000000000019
# ╟─a0000020-0000-0000-0000-000000000020
# ╠═a0000021-0000-0000-0000-000000000021
# ╟─a0000022-0000-0000-0000-000000000022
# ╠═a0000023-0000-0000-0000-000000000023
# ╟─a0000024-0000-0000-0000-000000000024
# ╠═a0000025-0000-0000-0000-000000000025
# ╟─a0000026-0000-0000-0000-000000000026
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
