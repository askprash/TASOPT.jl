### A Pluto.jl notebook ###
# v0.20.4
#
# TASOPT Interactive Comparison Dashboard (WGLMakie / Pluto edition)
#
# Run with:  julia --project=.. -e 'using Pluto; Pluto.run(notebook="interactive_dashboard_pluto.jl")'
#
# NOTE: The standalone GLMakie version (interactive_dashboard.jl) provides a
# richer native-window experience.  Use this notebook version when you need
# browser access or are on a headless server.

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity
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

# ╔═╡ a1000001-0000-0000-0000-000000000001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.add(["WGLMakie", "PlutoUI"])
    Pkg.instantiate()
end

# ╔═╡ a1000002-0000-0000-0000-000000000002
begin
    using WGLMakie
    using TASOPT
    using PlutoUI
    using Printf
    include(TASOPT.__TASOPTindices__)
    WGLMakie.activate!()
    TableOfContents(title="✈ TASOPT Comparison")
end

# ╔═╡ a1000003-0000-0000-0000-000000000003
md"""
# ✈ TASOPT Interactive Comparison Dashboard

> **Tip:** For the full interactive native experience (hover tooltips, click-to-switch, real-time sliders), run the standalone GLMakie script:
> ```julia
> julia --project=.. interactive_dashboard.jl ac1.jld2 ac2.jld2
> ```
> This Pluto version provides the same analyses with Pluto's slider-based interactivity.
"""

# ╔═╡ a1000004-0000-0000-0000-000000000004
md"""
## 📂 Load Aircraft

| | Path |
|:--|:--|
| **Aircraft 1** | $(@bind ac1_path TextField(80; default=joinpath(TASOPT.__TASOPTroot__, "IO/default_quicksave_aircraft.jld2"))) |
| **Aircraft 2** | $(@bind ac2_path TextField(80; default=joinpath(TASOPT.__TASOPTroot__, "IO/default_quicksave_aircraft.jld2"))) |

> Save a sized aircraft with `TASOPT.quicksave_aircraft(ac, "path.jld2")`
"""

# ╔═╡ a1000005-0000-0000-0000-000000000005
begin
    function try_load(path)
        isfile(path) || return (nothing, "File not found: $path")
        try
            ac = TASOPT.quickload_aircraft(path)
            ac.is_sized[1] || return (nothing, "Not sized")
            return (ac, nothing)
        catch e
            return (nothing, "Error: $e")
        end
    end
    ac1, err1 = try_load(ac1_path)
    ac2, err2 = try_load(ac2_path)
    nothing
end;

# ╔═╡ a1000006-0000-0000-0000-000000000006
if isnothing(ac1) || isnothing(ac2)
    md"""!!! danger "Not loaded"
        $(isnothing(ac1) ? "**AC1:** $err1" : "")
        $(isnothing(ac2) ? "**AC2:** $err2" : "")"""
else
    md"""!!! success "Both loaded"
        **AC1:** $(ac1.name)  |  **AC2:** $(ac2.name)"""
end

# ╔═╡ a1000007-0000-0000-0000-000000000007
md"""
## 🎛 Mission Point Selector

Use the slider below to select a mission point. All reactive plots update.

Mission point: $(@bind sel_pt_raw PlutoUI.Slider(1:iptotal, default=ipcruise1, show_value=true))
"""

# ╔═╡ a1000008-0000-0000-0000-000000000008
sel_pt = sel_pt_raw   # alias so we can reference below

# ╔═╡ a1000009-0000-0000-0000-000000000009
md"""Current mission point: **$(ip_labels[sel_pt])** (pt $(sel_pt))"""

# ╔═╡ a1000010-0000-0000-0000-000000000010
md"""
---
## ✏️ Aircraft Outlines (Top View)

Hover over regions to see tooltips.  In the native GLMakie app you can also click components.
"""

# ╔═╡ a1000011-0000-0000-0000-000000000011
# ── Colour helpers ──────────────────────────────────────────────────────────
begin
    C1_wgl = RGBf(0.27, 0.51, 0.71)
    C2_wgl = RGBf(0.85, 0.33, 0.10)
    AX_BG_wgl = RGBf(0.10, 0.10, 0.13)
    TXT_COL_wgl = RGBf(0.90, 0.90, 0.90)
end;

# ╔═╡ a1000012-0000-0000-0000-000000000012
# ── Wing polygon helper ─────────────────────────────────────────────────────
function wing_polygon_wgl(ac)
    w = ac.wing
    co = w.layout.root_chord
    λs = w.inboard.λ;  cs = co * λs
    λt = w.outboard.λ; ct = co * λt
    bo = w.layout.root_span; bs = w.layout.break_span; b = w.layout.span
    sw = w.layout.sweep; dx = w.layout.box_x
    tanL = tan(sw * π / 180.0); xax = 0.40
    xs = tanL*(bs-bo)/2; xt = tanL*(b-bo)/2
    xw = [co*(-xax)+dx, xs+cs*(-xax)+dx, xt+ct*(-xax)+dx,
          xt+ct*(1-xax)+dx, xs+cs*(1-xax)+dx, co*(1-xax)+dx]
    yw = [bo/2, bs/2, b/2, b/2, bs/2, bo/2]
    xs_all = vcat(xw, reverse(xw))
    ys_all = vcat(yw, reverse(-yw))
    return [Point2f(x, y) for (x, y) in zip(xs_all, ys_all)]
end;

# ╔═╡ a1000013-0000-0000-0000-000000000013
# ── Fuselage polygon helper ─────────────────────────────────────────────────
function fuselage_polygon_wgl(ac)
    f = ac.fuselage
    Rf=f.layout.radius; wfb=f.layout.bubble_center_y_offset
    an=f.layout.nose_radius; bt=f.layout.tail_radius
    x0=f.layout.x_nose; x1=f.layout.x_start_cylinder
    x2=f.layout.x_end_cylinder; xe=f.layout.x_end
    hw=Rf+wfb
    dy = startswith(string(f.layout.opt_tapers_to),"p") ? -hw : -0.2*hw
    N=20
    xf=zeros(N+N+1); yf=zeros(N+N+1)
    for i in 1:N
        t=(i-1)/(N-1); fx=cos(0.5π*t)
        xf[i]=x1+(x0-x1)*fx; yf[i]=hw*(1-fx^an)^(1/an)
    end
    for i in 1:N
        t=(i-1)/(N-1)
        xf[N+i]=x2+(xe-x2)*t; yf[N+i]=hw+dy*t^bt
    end
    xf[end]=xf[end-1]; yf[end]=0.0
    return [Point2f(x,y) for (x,y) in zip(vcat(xf,reverse(xf)), vcat(yf,reverse(-yf)))]
end;

# ╔═╡ a1000014-0000-0000-0000-000000000014
# ── HTail polygon helper ────────────────────────────────────────────────────
function htail_polygon_wgl(ac)
    h=ac.htail; f=ac.fuselage
    Sh=h.layout.S; ARh=h.layout.AR; λh=h.outboard.λ
    boh=h.layout.root_span; swh=h.layout.sweep; dx=h.layout.box_x
    bh=sqrt(Sh*ARh); coh=Sh/(boh+(bh-boh)*0.5*(1+λh)); cth=coh*λh
    tanLh=tan(swh*π/180); xax=0.40
    xoLE=coh*(-xax)+dx; xoTE=coh*(1-xax)+dx
    xtLE=cth*(-xax)+dx+0.5*(bh-boh)*tanLh; xtTE=cth*(1-xax)+dx+0.5*(bh-boh)*tanLh
    if startswith(string(f.layout.opt_tapers_to),"p")
        xcLE=xoLE; xcTE=xoTE; ycLE=0.5*boh; ycTE=0.5*boh
    else
        xcLE=coh*(-xax)+dx-0.5*boh*tanLh; xcTE=coh*(1-xax)+dx-0.5*boh*tanLh
        ycLE=0.0; ycTE=0.0
    end
    xh=[xcLE,xoLE,xtLE,xtTE,xoTE,xcTE]; yh=[ycLE,0.5*boh,0.5*bh,0.5*bh,0.5*boh,ycTE]
    return [Point2f(x,y) for (x,y) in zip(vcat(xh,reverse(xh)), vcat(yh,reverse(-yh)))]
end;

# ╔═╡ a1000015-0000-0000-0000-000000000015
if !isnothing(ac1) && !isnothing(ac2)
    fig_outline = Figure(size=(1200, 500), backgroundcolor=RGBf(0.12,0.12,0.15))

    for (col_idx, (ac, col)) in enumerate([(ac1, C1_wgl), (ac2, C2_wgl)])
        ax = Axis(fig_outline[1, col_idx];
            title=ac.name, titlecolor=col, titlefontsize=13,
            xlabel="x [m]", ylabel="y [m]",
            aspect=DataAspect(), backgroundcolor=AX_BG_wgl,
            xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
            xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)

        # Wing
        poly!(ax, wing_polygon_wgl(ac); color=(col, 0.18), strokecolor=col, strokewidth=2.5,
            inspector_label=(_, _, _) -> begin
                w = ac.wing
                @sprintf("Wing: span=%.1fm AR=%.2f sweep=%.1f° S=%.0fm²\nCruise CL=%.4f  CD_wing=%.5f",
                    w.layout.span, w.layout.AR, w.layout.sweep, w.layout.S,
                    ac.para[iaCL,ipcruise1,1], ac.para[iaCDwing,ipcruise1,1])
            end)

        # Fuselage
        poly!(ax, fuselage_polygon_wgl(ac);
            color=(RGBf(0.7,0.7,0.7), 0.55), strokecolor=col, strokewidth=2.5,
            inspector_label=(_, _, _) ->
                @sprintf("Fuselage: L=%.1fm R=%.3fm W=%.2ft",
                    ac.fuselage.layout.x_end-ac.fuselage.layout.x_nose,
                    ac.fuselage.layout.radius, ac.fuselage.weight/9.81e3))

        # HTail
        poly!(ax, htail_polygon_wgl(ac);
            color=(col, 0.30), strokecolor=col, strokewidth=2.0,
            inspector_label=(_, _, _) ->
                @sprintf("HTail: S=%.1fm² AR=%.2f sweep=%.1f° W=%.3ft",
                    ac.htail.layout.S, ac.htail.layout.AR,
                    ac.htail.layout.sweep, ac.htail.weight/9.81e3))

        # Engines
        pg=ac.parg; w=ac.wing; D=pg[igdfan]; lnace=pg[iglnace]
        b=w.layout.span; bo=w.layout.root_span
        ηs=w.layout.ηs; tanL=tan(w.layout.sweep*π/180); dx=w.layout.box_x
        neng=Int(pg[igneng])
        yi = neng==2 ? [ηs*b/2] : collect(range(bo/2+2D, b/2*3/4; length=Int(neng/2)))
        for y in yi
            for sg in (-1,1)
                λs=w.inboard.λ; λt=w.outboard.λ; co=w.layout.root_chord
                η=y/(b/2); ηo=bo/b
                ci = η<=ηs ? co*(1+(λs-1)*(η-ηo)/(ηs-ηo)) : co*(λs+(λt-λs)*(η-ηs)/(1-ηs))
                xi=tanL*(y-bo/2)-0.4*ci+dx-1.0; yc=sg*y
                rect=[Point2f(xi,yc-D/2),Point2f(xi+lnace,yc-D/2),
                      Point2f(xi+lnace,yc+D/2),Point2f(xi,yc+D/2)]
                poly!(ax, rect; color=(RGBf(0.9,0.2,0.1),0.35),
                    strokecolor=RGBf(0.9,0.2,0.1), strokewidth=1.5,
                    inspector_label=(_, _, _) ->
                        @sprintf("Engine: BPR=%.2f OPR=%.1f Tt4=%.0fK TSFC=%.4fmg/N/s",
                            ac.pare[ieBPR,ipcruise1,1], ac.pare[ieOPR,ipcruise1,1],
                            ac.pare[ieTt4,ipcruise1,1], ac.pare[ieTSFC,ipcruise1,1]*1e6))
            end
        end

        # CG & NP
        xCGf=pg[igxCGfwd]; xCGa=pg[igxCGaft]; xNP=pg[igxNP]
        lines!(ax,[xCGf,xCGa],[0.0,0.0]; color=:white, linewidth=3)
        scatter!(ax,[xNP],[0.0]; color=RGBf(1.0,0.85,0.0), marker=:utriangle, markersize=12)
        text!(ax, xNP, 1.5; text="NP", color=RGBf(1,0.85,0), fontsize=9, align=(:center,:bottom))
    end

    DataInspector(fig_outline;
        indicator_color=RGBf(1,0.85,0), fontsize=11,
        background_color=RGBAf(0.1,0.1,0.13,0.92),
        text_color=TXT_COL_wgl)

    fig_outline
end

# ╔═╡ a1000016-0000-0000-0000-000000000016
md"""
---
## 📈 Mission Profile

All panels update when you move the **Mission Point** slider at the top.
"""

# ╔═╡ a1000017-0000-0000-0000-000000000017
if !isnothing(ac1) && !isnothing(ac2)
    fig_mission = Figure(size=(1100, 580), backgroundcolor=RGBf(0.12,0.12,0.15))
    pts = collect(1:iptotal)

    for (row, ylabel, data1, data2, title) in [
        (1, "Alt [km]",  ac1.para[iaalt,:,1]./1e3,          ac2.para[iaalt,:,1]./1e3,          "Altitude"),
        (1, "Mach",      ac1.para[iaMach,:,1],               ac2.para[iaMach,:,1],               "Mach"),
        (2, "L/D",       ac1.para[iaCL,:,1]./ac1.para[iaCD,:,1], ac2.para[iaCL,:,1]./ac2.para[iaCD,:,1], "L/D"),
        (2, "ROC [m/s]", ac1.para[iaROC,:,1],                ac2.para[iaROC,:,1],                "Rate of Climb"),
    ]
        col = (row == 1 && ylabel == "Alt [km]") ? 1 :
              (row == 1 && ylabel == "Mach") ? 2 :
              (row == 2 && ylabel == "L/D") ? 1 : 2
        ax = Axis(fig_mission[row, col];
            title=title, ylabel=ylabel, xlabel=(row==2 ? "Mission Pt" : ""),
            backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
            xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
            xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)
        lines!(ax, pts, data1; color=C1_wgl, linewidth=2.5, label=ac1.name)
        lines!(ax, pts, data2; color=C2_wgl, linewidth=2.5, label=ac2.name)
        vlines!(ax, [sel_pt]; color=(RGBf(1,0.85,0),0.9), linewidth=2.5, linestyle=:dash)
        row == 1 && col == 1 && axislegend(ax; position=:rb, labelsize=9)
    end
    fig_mission
end

# ╔═╡ a1000018-0000-0000-0000-000000000018
md"""
---
## ⚙ Engine Analysis

**Station diagram** at the selected mission point + trends along the mission.
"""

# ╔═╡ a1000019-0000-0000-0000-000000000019
if !isnothing(ac1) && !isnothing(ac2)
    Tvars = [ieTt0,ieTt2,ieTt21,ieTt25,ieTt3,ieTt4,ieTt41,ieTt45,ieTt5,ieTt7]
    Tnames= ["Tt₀","Tt₂","Tt₂₁","Tt₂₅","Tt₃","Tt₄","Tt₄₁","Tt₄₅","Tt₅","Tt₇"]
    Pvars = [iept0,iept2,iept21,iept25,iept3,iept4,iept41,iept45,iept5,iept7]
    Pnames= ["pt₀","pt₂","pt₂₁","pt₂₅","pt₃","pt₄","pt₄₁","pt₄₅","pt₅","pt₇"]
    ns    = 1:length(Tvars)

    T1 = [ac1.pare[v, sel_pt, 1] for v in Tvars]
    T2 = [ac2.pare[v, sel_pt, 1] for v in Tvars]
    P1 = [ac1.pare[v, sel_pt, 1]/1e3 for v in Pvars]
    P2 = [ac2.pare[v, sel_pt, 1]/1e3 for v in Pvars]

    fig_eng = Figure(size=(1100, 700), backgroundcolor=RGBf(0.12,0.12,0.15))

    ax_T = Axis(fig_eng[1, 1:2];
        title="Engine Station Temperatures – $(ip_labels[sel_pt])",
        xticks=(collect(ns), Tnames), ylabel="Tt [K]",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        ylabelcolor=TXT_COL_wgl)
    lines!(ax_T, collect(ns), T1; color=C1_wgl, linewidth=2.5, marker=:circle, markersize=9, label=ac1.name)
    lines!(ax_T, collect(ns), T2; color=C2_wgl, linewidth=2.5, marker=:circle, markersize=9, label=ac2.name)
    axislegend(ax_T; position=:lt, labelsize=10)

    ax_P = Axis(fig_eng[2, 1:2];
        title="Engine Station Total Pressures",
        xticks=(collect(1:length(Pvars)), Pnames), ylabel="pt [kPa]",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        ylabelcolor=TXT_COL_wgl)
    lines!(ax_P, collect(1:length(Pvars)), P1; color=C1_wgl, linewidth=2.5, marker=:circle, markersize=9)
    lines!(ax_P, collect(1:length(Pvars)), P2; color=C2_wgl, linewidth=2.5, marker=:circle, markersize=9)

    # Mission trends
    pts = collect(1:iptotal)
    for (col_idx, (ylabel, d1, d2, ttl)) in enumerate([
        ("Tt4 [K]",       ac1.pare[ieTt4,:,1],          ac2.pare[ieTt4,:,1],       "Tt4 Along Mission"),
        ("TSFC [mg/N/s]", ac1.pare[ieTSFC,:,1].*1e6,    ac2.pare[ieTSFC,:,1].*1e6, "TSFC Along Mission"),
        ("BPR",           ac1.pare[ieBPR,:,1],           ac2.pare[ieBPR,:,1],       "BPR Along Mission"),
        ("OPR",           ac1.pare[ieOPR,:,1],           ac2.pare[ieOPR,:,1],       "OPR Along Mission"),
    ])
        r = col_idx <= 2 ? 3 : 4
        c = mod(col_idx-1, 2) + 1
        ax = Axis(fig_eng[r, c];
            title=ttl, ylabel=ylabel, xlabel="Mission Pt",
            backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
            xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
            ylabelcolor=TXT_COL_wgl, xlabelcolor=TXT_COL_wgl)
        lines!(ax, pts, d1; color=C1_wgl, linewidth=2.5)
        lines!(ax, pts, d2; color=C2_wgl, linewidth=2.5)
        vlines!(ax, [sel_pt]; color=(RGBf(1,0.85,0),0.9), linewidth=2.5, linestyle=:dash)
    end

    fig_eng
end

# ╔═╡ a1000020-0000-0000-0000-000000000020
md"""
---
## 💨 Drag Breakdown

Component drag breakdown **at the selected mission point** + trends along the mission.
"""

# ╔═╡ a1000021-0000-0000-0000-000000000021
if !isnothing(ac1) && !isnothing(ac2)
    drag_keys = ["CDi","CDfuse","CDwing","CDhtail","CDvtail","CDnace"]
    drag_inds = [iaCDi, iaCDfuse, iaCDwing, iaCDhtail, iaCDvtail, iaCDnace]
    drag_cols_wgl = [RGBf(0.27,0.51,0.71), RGBf(0.85,0.33,0.10), RGBf(0.18,0.63,0.22),
                     RGBf(0.95,0.60,0.07), RGBf(0.60,0.20,0.70), RGBf(0.90,0.40,0.55)]
    S1 = ac1.wing.layout.S; S2 = ac2.wing.layout.S

    cd1_pt = [ac1.para[idx, sel_pt, 1]*S1*1e4 for idx in drag_inds]
    cd2_pt = [ac2.para[idx, sel_pt, 1]*S2*1e4 for idx in drag_inds]
    nx = length(drag_keys)

    fig_drag = Figure(size=(1100, 650), backgroundcolor=RGBf(0.12,0.12,0.15))
    ax_bar = Axis(fig_drag[1, 1:2];
        title="Drag Area (CD×S) at $(ip_labels[sel_pt])  [cm²]",
        xticks=(1:nx, drag_keys), ylabel="CD×S [cm²]",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        ylabelcolor=TXT_COL_wgl)
    barplot!(ax_bar, (1:nx) .- 0.22, cd1_pt;
        width=0.38, color=[(c, 0.85) for c in drag_cols_wgl], label=ac1.name)
    barplot!(ax_bar, (1:nx) .+ 0.22, cd2_pt;
        width=0.38, color=[(c, 0.55) for c in drag_cols_wgl],
        strokecolor=drag_cols_wgl, strokewidth=1.0, label=ac2.name)
    axislegend(ax_bar; position=:rt, labelsize=10)

    # L/D and CD trends
    pts = collect(1:iptotal)
    CD1all=ac1.para[iaCD,:,1]; CD2all=ac2.para[iaCD,:,1]
    LoD1=ac1.para[iaCL,:,1]./CD1all; LoD2=ac2.para[iaCL,:,1]./CD2all

    ax_LoD = Axis(fig_drag[2,1]; title="L/D", ylabel="L/D", xlabel="Mission Pt",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)
    lines!(ax_LoD, pts, LoD1; color=C1_wgl, linewidth=2.5, label=ac1.name)
    lines!(ax_LoD, pts, LoD2; color=C2_wgl, linewidth=2.5, label=ac2.name)
    vlines!(ax_LoD, [sel_pt]; color=(RGBf(1,0.85,0),0.9), linewidth=2.5, linestyle=:dash)
    axislegend(ax_LoD; position=:lb, labelsize=9)

    ax_CD = Axis(fig_drag[2,2]; title="Total CD [×10⁻⁴]", ylabel="CD×10⁴", xlabel="Mission Pt",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)
    lines!(ax_CD, pts, CD1all.*1e4; color=C1_wgl, linewidth=2.5)
    lines!(ax_CD, pts, CD2all.*1e4; color=C2_wgl, linewidth=2.5)
    vlines!(ax_CD, [sel_pt]; color=(RGBf(1,0.85,0),0.9), linewidth=2.5, linestyle=:dash)

    fig_drag
end

# ╔═╡ a1000022-0000-0000-0000-000000000022
md"""
---
## ⚖ Weight Breakdown
"""

# ╔═╡ a1000023-0000-0000-0000-000000000023
if !isnothing(ac1) && !isnothing(ac2)
    function get_weights_wgl(ac)
        pg=ac.parg; w=ac.wing; lg=ac.landing_gear
        Wbox=w.inboard.webs.weight.W+w.inboard.caps.weight.W
        Wwing=Wbox*(1+w.weight_frac_flap+w.weight_frac_slat+w.weight_frac_ailerons+
                    w.weight_frac_leading_trailing_edge+w.weight_frac_ribs+
                    w.weight_frac_spoilers+w.weight_frac_attachments)
        Wadd=pg[igWMTO]*ac.fuselage.HPE_sys.W+lg.nose_gear.weight.W+lg.main_gear.weight.W
        lbls=["Fuselage","Wing","H-Tail","V-Tail","Engine Sys","Fuel Tank","Misc"]
        vals=[ac.fuselage.weight,Wwing,ac.htail.weight,ac.vtail.weight,
              pg[igWtesys],pg[igWftank],Wadd]./9.81e3
        return lbls, vals
    end
    wlabels,wv1=get_weights_wgl(ac1); _,wv2=get_weights_wgl(ac2)
    nx=length(wlabels)
    w_cols=[RGBf(0.27,0.51,0.71),RGBf(0.18,0.63,0.22),RGBf(0.95,0.60,0.07),
            RGBf(0.60,0.20,0.70),RGBf(0.85,0.33,0.10),RGBf(0.90,0.40,0.55),RGBf(0.50,0.50,0.55)]

    fig_w = Figure(size=(1000, 400), backgroundcolor=RGBf(0.12,0.12,0.15))
    ax_w = Axis(fig_w[1,1];
        title="Structural Weight Breakdown [t]",
        xticks=(1:nx, wlabels), ylabel="Weight [t]",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        ylabelcolor=TXT_COL_wgl)
    barplot!(ax_w, (1:nx).-0.22, wv1; width=0.38,
        color=[(c,0.85) for c in w_cols], label=ac1.name)
    barplot!(ax_w, (1:nx).+0.22, wv2; width=0.38,
        color=[(c,0.55) for c in w_cols], strokecolor=w_cols, strokewidth=1.0, label=ac2.name)
    axislegend(ax_w; position=:rt, labelsize=10)
    fig_w
end

# ╔═╡ a1000024-0000-0000-0000-000000000024
md"""
---
## 📐 Wing & Geometry
"""

# ╔═╡ a1000025-0000-0000-0000-000000000025
if !isnothing(ac1) && !isnothing(ac2)
    function chord_dist_wgl(ac)
        w=ac.wing; bo=w.layout.root_span; bs=w.layout.break_span; b=w.layout.span
        co=w.layout.root_chord; cs=co*w.inboard.λ; ct=co*w.outboard.λ
        [0.0,bo/b,bs/b,1.0], [co,co,cs,ct]
    end
    η1,c1=chord_dist_wgl(ac1); η2,c2=chord_dist_wgl(ac2)

    fig_wing = Figure(size=(1100, 520), backgroundcolor=RGBf(0.12,0.12,0.15))

    ax_c = Axis(fig_wing[1,1]; title="Wing Chord Distribution",
        xlabel="η=2y/b", ylabel="Chord [m]",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)
    lines!(ax_c, η1, c1; color=C1_wgl, linewidth=2.5, marker=:circle, markersize=9, label=ac1.name)
    lines!(ax_c, η2, c2; color=C2_wgl, linewidth=2.5, marker=:circle, markersize=9, label=ac2.name)
    axislegend(ax_c; position=:tr, labelsize=10)

    ax_f = Axis(fig_wing[1,2]; title="Fuselage Cross-Section",
        xlabel="y [m]", ylabel="z [m]", aspect=DataAspect(),
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        xlabelcolor=TXT_COL_wgl, ylabelcolor=TXT_COL_wgl)
    R1=ac1.fuselage.layout.radius; R2=ac2.fuselage.layout.radius
    θ=LinRange(0,2π,200)
    lines!(ax_f, R1.*cos.(θ), R1.*sin.(θ); color=C1_wgl, linewidth=2.5,
        label="$(ac1.name) R=$(round(R1,digits=3))m")
    lines!(ax_f, R2.*cos.(θ), R2.*sin.(θ); color=C2_wgl, linewidth=2.5,
        linestyle=:dash, label="$(ac2.name) R=$(round(R2,digits=3))m")
    axislegend(ax_f; position=:lb, labelsize=9)

    Lf1=ac1.fuselage.layout.x_end-ac1.fuselage.layout.x_nose
    Lf2=ac2.fuselage.layout.x_end-ac2.fuselage.layout.x_nose
    glabels=["Span [m]","Fuse L [m]","S [m²]","Root c [m]","AR","Sweep [°]"]
    gv1=[ac1.wing.layout.span,Lf1,ac1.wing.layout.S,
         ac1.wing.layout.root_chord,ac1.wing.layout.AR,ac1.wing.layout.sweep]
    gv2=[ac2.wing.layout.span,Lf2,ac2.wing.layout.S,
         ac2.wing.layout.root_chord,ac2.wing.layout.AR,ac2.wing.layout.sweep]
    gnx=length(glabels)

    ax_g = Axis(fig_wing[2,1:2]; title="Key Geometry Comparison",
        xticks=(1:gnx, glabels), ylabel="Value",
        backgroundcolor=AX_BG_wgl, titlecolor=TXT_COL_wgl,
        xticklabelcolor=TXT_COL_wgl, yticklabelcolor=TXT_COL_wgl,
        ylabelcolor=TXT_COL_wgl)
    barplot!(ax_g, (1:gnx).-0.22, gv1; width=0.38, color=(C1_wgl,0.85), label=ac1.name)
    barplot!(ax_g, (1:gnx).+0.22, gv2; width=0.38, color=(C2_wgl,0.85), label=ac2.name)
    axislegend(ax_g; position=:tr, labelsize=10)

    fig_wing
end

# ╔═╡ a1000026-0000-0000-0000-000000000026
md"""
---
*TASOPT Interactive Comparison Dashboard — WGLMakie/Pluto edition*

For the full interactive native experience with hover tooltips and click events, run:
```julia
julia --project=.. interactive_dashboard.jl ac1.jld2 ac2.jld2
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
PlutoUI = "0.7"
WGLMakie = "0.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = ""

# ╔═╡ Cell Order:
# ╠═a1000001-0000-0000-0000-000000000001
# ╟─a1000002-0000-0000-0000-000000000002
# ╟─a1000003-0000-0000-0000-000000000003
# ╟─a1000004-0000-0000-0000-000000000004
# ╠═a1000005-0000-0000-0000-000000000005
# ╟─a1000006-0000-0000-0000-000000000006
# ╟─a1000007-0000-0000-0000-000000000007
# ╠═a1000008-0000-0000-0000-000000000008
# ╟─a1000009-0000-0000-0000-000000000009
# ╟─a1000010-0000-0000-0000-000000000010
# ╠═a1000011-0000-0000-0000-000000000011
# ╠═a1000012-0000-0000-0000-000000000012
# ╠═a1000013-0000-0000-0000-000000000013
# ╠═a1000014-0000-0000-0000-000000000014
# ╠═a1000015-0000-0000-0000-000000000015
# ╟─a1000016-0000-0000-0000-000000000016
# ╠═a1000017-0000-0000-0000-000000000017
# ╟─a1000018-0000-0000-0000-000000000018
# ╠═a1000019-0000-0000-0000-000000000019
# ╟─a1000020-0000-0000-0000-000000000020
# ╠═a1000021-0000-0000-0000-000000000021
# ╟─a1000022-0000-0000-0000-000000000022
# ╠═a1000023-0000-0000-0000-000000000023
# ╟─a1000024-0000-0000-0000-000000000024
# ╠═a1000025-0000-0000-0000-000000000025
# ╟─a1000026-0000-0000-0000-000000000026
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
