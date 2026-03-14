##############################################################################
#  TASOPT Interactive Comparison Dashboard
#  interactive_dashboard.jl
#
#  Launches a native GLMakie window that lets you interactively compare two
#  sized TASOPT aircraft side-by-side.
#
#  FEATURES
#  --------
#  • Aircraft top-view outlines with hover-tooltips (via DataInspector)
#      – Wing    → span, AR, sweep, chord, cruise CL/CD
#      – Fuselage → length, radius, weight
#      – H-Tail  → span, area, sweep, weight
#      – Engines → BPR, OPR, fan diameter, Tt4, TSFC at cruise
#  • Click a component → jump to the relevant analysis tab
#  • Mission-point slider → all analyses update live
#  • Five analysis tabs (rebuilt on switch):
#      1. Mission Profile  – altitude, L/D, ROC, weight fraction along mission
#      2. Engine Analysis  – station Tt/pt diagram + Tt4/TSFC/BPR trends
#      3. Drag Breakdown   – component drag bars + CD / L/D trends vs mission
#      4. Weights          – structural breakdown + MTOW pie charts
#      5. Wing & Geometry  – chord distribution, fuselage cross-section, key dims
#
#  USAGE
#  -----
#  From the shell (with a display / DISPLAY set):
#
#      julia --project=.. interactive_dashboard.jl ac1.jld2 ac2.jld2
#
#  From the Julia REPL:
#
#      include("interactive_dashboard.jl")
#      launch_dashboard(ac1, ac2)          # ac1/ac2 are TASOPT aircraft objects
#
#  Save an aircraft with:  TASOPT.quicksave_aircraft(ac, "path/to/file.jld2")
##############################################################################

using GLMakie
using TASOPT
using Printf

include(TASOPT.__TASOPTindices__)

# ─────────────────────────────────────────────────────────────────────────────
# Theme / colours
# ─────────────────────────────────────────────────────────────────────────────
const C1       = RGBf(0.27, 0.51, 0.71)   # steelblue  – aircraft 1
const C2       = RGBf(0.85, 0.33, 0.10)   # tomato     – aircraft 2
const DARK_BG  = RGBf(0.12, 0.12, 0.15)
const MID_BG   = RGBf(0.18, 0.18, 0.22)
const AX_BG    = RGBf(0.10, 0.10, 0.13)
const TXT_COL  = RGBf(0.90, 0.90, 0.90)
const YELLOW   = RGBf(1.0,  0.85, 0.0)

set_theme!(
    backgroundcolor = DARK_BG,
    textcolor       = TXT_COL,
    Axis = (
        backgroundcolor = AX_BG,
        xgridcolor = RGBAf(1,1,1,0.06),
        ygridcolor = RGBAf(1,1,1,0.06),
        xticklabelcolor = TXT_COL,
        yticklabelcolor = TXT_COL,
        xlabelcolor     = TXT_COL,
        ylabelcolor     = TXT_COL,
        titlecolor      = TXT_COL,
        leftspinecolor  = RGBAf(1,1,1,0.3),
        bottomspinecolor= RGBAf(1,1,1,0.3),
        rightspinecolor = :transparent,
        topspinecolor   = :transparent,
    ),
    Legend = (
        backgroundcolor = RGBAf(0.1,0.1,0.13,0.85),
        framecolor      = RGBAf(1,1,1,0.2),
        labelcolor      = TXT_COL,
    ),
)

# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers – replicate stickfig logic but return Point2f arrays
# ─────────────────────────────────────────────────────────────────────────────

"""Return closed polygon for the FULL wing (upper + lower half, top view)."""
function wing_polygon(ac)
    w   = ac.wing
    co  = w.layout.root_chord
    λs  = w.inboard.λ;   cs = co * λs
    λt  = w.outboard.λ;  ct = co * λt
    bo  = w.layout.root_span
    bs  = w.layout.break_span
    b   = w.layout.span
    sw  = w.layout.sweep
    dx  = w.layout.box_x
    tanL = tan(sw * π / 180.0)
    xax  = 0.40

    xs  = tanL * (bs - bo) / 2.0
    xt  = tanL * (b  - bo) / 2.0

    # Upper half: root-LE → break-LE → tip-LE → tip-TE → break-TE → root-TE
    xw = [co * (-xax) + dx,
          xs + cs * (-xax) + dx,
          xt + ct * (-xax) + dx,
          xt + ct * (1-xax) + dx,
          xs + cs * (1-xax) + dx,
             co * (1-xax) + dx]
    yw = [bo/2, bs/2, b/2, b/2, bs/2, bo/2]

    # Full (upper ++ reversed lower)
    xs_all = vcat(xw, reverse(xw))
    ys_all = vcat(yw, reverse(-yw))
    return [Point2f(x, y) for (x, y) in zip(xs_all, ys_all)]
end

"""Return closed polygon for the fuselage outline (top view)."""
function fuselage_polygon(ac)
    f  = ac.fuselage
    Rf = f.layout.radius
    wfb= f.layout.bubble_center_y_offset
    an = f.layout.nose_radius
    bt = f.layout.tail_radius
    x0 = f.layout.x_nose
    x1 = f.layout.x_start_cylinder
    x2 = f.layout.x_end_cylinder
    xe = f.layout.x_end
    hw = Rf + wfb
    dy = startswith(string(f.layout.opt_tapers_to), "p") ? -hw : -0.2*hw
    N  = 20

    xf = zeros(N + N + 1)
    yf = zeros(N + N + 1)
    for i in 1:N
        t  = (i-1)/(N-1); fx = cos(0.5π*t)
        xf[i] = x1 + (x0 - x1)*fx
        yf[i] = hw * (1 - fx^an)^(1/an)
    end
    for i in 1:N
        t = (i-1)/(N-1)
        xf[N+i] = x2 + (xe - x2)*t
        yf[N+i] = hw + dy*t^bt
    end
    xf[end] = xf[end-1]; yf[end] = 0.0

    xs = vcat(xf, reverse(xf))
    ys = vcat(yf, reverse(-yf))
    return [Point2f(x, y) for (x, y) in zip(xs, ys)]
end

"""Return closed polygon for the horizontal tail (top view)."""
function htail_polygon(ac)
    h  = ac.htail
    f  = ac.fuselage
    Sh = h.layout.S; ARh = h.layout.AR; λh = h.outboard.λ
    boh= h.layout.root_span; swh= h.layout.sweep; dx= h.layout.box_x
    bh = sqrt(Sh * ARh)
    coh= Sh / (boh + (bh - boh)*0.5*(1+λh))
    cth= coh * λh
    tanLh = tan(swh * π/180.0)
    xax  = 0.40

    xoLE = coh*(-xax) + dx
    xoTE = coh*(1-xax) + dx
    xtLE = cth*(-xax) + dx + 0.5*(bh-boh)*tanLh
    xtTE = cth*(1-xax) + dx + 0.5*(bh-boh)*tanLh

    if startswith(string(f.layout.opt_tapers_to), "p")
        xcLE = xoLE; xcTE = xoTE; ycLE = 0.5*boh; ycTE = 0.5*boh
    else
        xcLE = coh*(-xax) + dx - 0.5*boh*tanLh
        xcTE = coh*(1-xax) + dx - 0.5*boh*tanLh
        ycLE = 0.0; ycTE = 0.0
    end

    xh = [xcLE, xoLE, xtLE, xtTE, xoTE, xcTE]
    yh = [ycLE, 0.5*boh, 0.5*bh, 0.5*bh, 0.5*boh, ycTE]

    xs = vcat(xh, reverse(xh))
    ys = vcat(yh, reverse(-yh))
    return [Point2f(x, y) for (x, y) in zip(xs, ys)]
end

"""Return vector of nacelle rectangles (one per engine, both sides)."""
function engine_rects(ac)
    pg = ac.parg; w = ac.wing
    D     = pg[igdfan]
    neng  = Int(pg[igneng])
    lnace = pg[iglnace]
    b = w.layout.span; bo = w.layout.root_span
    co= w.layout.root_chord; λs=w.inboard.λ; λt=w.outboard.λ
    ηs   = w.layout.ηs
    tanL = tan(w.layout.sweep * π/180.0)
    dx   = w.layout.box_x

    dy = 2D
    yi = (neng == 2) ? [ηs*b/2] :
                        collect(range(bo/2 + dy, b/2 * 3/4; length=Int(neng/2)))

    rects = Vector{Vector{Point2f}}()
    for y in yi
        η  = y / (b/2); ηo = bo/b
        ci = η <= ηs ? co*(1+(λs-1)*(η-ηo)/(ηs-ηo)) : co*(λs+(λt-λs)*(η-ηs)/(1-ηs))
        xi = tanL*(y - bo/2) - 0.4*ci + dx - 1.0
        for sg in (-1, 1)
            yc = sg * y
            push!(rects, [
                Point2f(xi,        yc - D/2),
                Point2f(xi+lnace,  yc - D/2),
                Point2f(xi+lnace,  yc + D/2),
                Point2f(xi,        yc + D/2),
            ])
        end
    end
    return rects
end

# ─────────────────────────────────────────────────────────────────────────────
# Draw a single aircraft on an Axis; return plot handles for pick()
# ─────────────────────────────────────────────────────────────────────────────

function draw_aircraft!(ax, ac, col)
    pg = ac.parg; para = view(ac.para, :, :, 1); pare = view(ac.pare, :, :, 1)

    # Wing ────────────────────────────────────────────────────────────────────
    w = ac.wing
    wp = wing_polygon(ac)
    wing_plt = poly!(ax, wp;
        color      = (col, 0.18),
        strokecolor = col,
        strokewidth = 2.5,
        inspector_label = (_, _, _) ->
            @sprintf("Wing\nSpan:     %.1f m     AR: %.2f\nSweep:    %.1f°     S:  %.0f m²\nRoot c:   %.2f m   Tip c: %.2f m\nλ_in: %.3f  λ_out: %.3f\n—— Cruise ——\nCL: %.4f   CD_wing: %.5f\nL/D: %.2f",
                w.layout.span, w.layout.AR, w.layout.sweep, w.layout.S,
                w.layout.root_chord, w.layout.root_chord * w.outboard.λ,
                w.inboard.λ, w.outboard.λ,
                para[iaCL, ipcruise1], para[iaCDwing, ipcruise1],
                para[iaCL, ipcruise1] / para[iaCD, ipcruise1])
    )

    # Panel break dashes (inner / outer)
    bo = w.layout.root_span; bs = w.layout.break_span; b = w.layout.span
    sw = w.layout.sweep; tanL = tan(sw*π/180); dx = w.layout.box_x
    xax = 0.40; co = w.layout.root_chord; λs = w.inboard.λ; λt = w.outboard.λ
    cs = co*λs; ct = co*λt
    xs = tanL*(bs-bo)/2; xt = tanL*(b-bo)/2
    for (x1, x2, yc) in [
        (xs + cs*(-xax)+dx, xs+cs*(1-xax)+dx, bs/2),
        (bo/2*(0) + dx - 0.4*co, dx + 0.6*co, bo/2),
    ]
        lines!(ax, [x1, x2], [yc, yc];  color=(col,0.35), linewidth=1, linestyle=:dash)
        lines!(ax, [x1, x2], [-yc,-yc]; color=(col,0.35), linewidth=1, linestyle=:dash)
    end

    # Fuselage ────────────────────────────────────────────────────────────────
    fp = fuselage_polygon(ac)
    fuse_plt = poly!(ax, fp;
        color       = (RGBf(0.7,0.7,0.7), 0.55),
        strokecolor  = col,
        strokewidth  = 2.5,
        inspector_label = (_, _, _) ->
            @sprintf("Fuselage\nLength:   %.1f m\nRadius:   %.3f m\nWeight:   %.2f t\nPax shell: %.1f - %.1f m\nnDecks: %d",
                ac.fuselage.layout.x_end - ac.fuselage.layout.x_nose,
                ac.fuselage.layout.radius,
                ac.fuselage.weight / 9.81e3,
                ac.fuselage.layout.x_pressure_shell_fwd,
                ac.fuselage.layout.x_pressure_shell_aft,
                ac.fuselage.n_decks)
    )

    # H-Tail ──────────────────────────────────────────────────────────────────
    hp = htail_polygon(ac)
    htail_plt = poly!(ax, hp;
        color       = (col, 0.30),
        strokecolor  = col,
        strokewidth  = 2.0,
        inspector_label = (_, _, _) ->
            @sprintf("H-Tail\nSpan: %.1f m   S: %.1f m2\nSweep: %.1f deg   AR: %.2f\nlambda_tip: %.3f\nWeight: %.3f t\nCL_h (cruise): %.4f",
                sqrt(ac.htail.layout.S * ac.htail.layout.AR),
                ac.htail.layout.S,
                ac.htail.layout.sweep, ac.htail.layout.AR,
                ac.htail.outboard.λ,
                ac.htail.weight / 9.81e3,
                para[iaCLh, ipcruise1])
    )

    # Engines ─────────────────────────────────────────────────────────────────
    eng_plts = map(engine_rects(ac)) do rect
        poly!(ax, rect;
            color       = (RGBf(0.9,0.2,0.1), 0.35),
            strokecolor  = RGBf(0.9,0.2,0.1),
            strokewidth  = 1.5,
            inspector_label = (_, _, _) ->
                @sprintf("Engine\nFan Dia:  %.3f m\nN_eng:    %d\n-- Cruise --\nBPR:  %.2f    OPR:  %.1f\nTt4:  %.0f K\nTSFC: %.4f mg/N/s\nFPR:  %.3f    LPC: %.3f   HPC: %.2f",
                    pg[igdfan], Int(pg[igneng]),
                    pare[ieBPR,  ipcruise1], pare[ieOPR,  ipcruise1],
                    pare[ieTt4,  ipcruise1],
                    pare[ieTSFC, ipcruise1]*1e6,
                    pare[iepif,  ipcruise1], pare[iepilc, ipcruise1], pare[iepihc, ipcruise1])
        )
    end

    # CG range & NP ───────────────────────────────────────────────────────────
    xCGf = pg[igxCGfwd]; xCGa = pg[igxCGaft]; xNP = pg[igxNP]
    lines!(ax, [xCGf, xCGa], [0.0, 0.0]; color=:white, linewidth=3)
    scatter!(ax, [xCGf, xCGa], [0.0, 0.0]; color=:white, marker=:vline, markersize=14)
    scatter!(ax, [xNP], [0.0]; color=YELLOW, marker=:utriangle, markersize=12)
    text!(ax, xNP, 1.5; text="NP", color=YELLOW, fontsize=9, align=(:center,:bottom))

    return (wing=wing_plt, fuselage=fuse_plt, htail=htail_plt, engines=eng_plts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Content panels (one per tab)
# ─────────────────────────────────────────────────────────────────────────────

function build_mission_panel!(gl, ac1, ac2, sel_pt)
    pts = collect(1:iptotal)

    alt1  = ac1.para[iaalt,  :, 1] ./ 1e3
    alt2  = ac2.para[iaalt,  :, 1] ./ 1e3
    LoD1  = ac1.para[iaCL, :, 1] ./ ac1.para[iaCD, :, 1]
    LoD2  = ac2.para[iaCL, :, 1] ./ ac2.para[iaCD, :, 1]
    roc1  = ac1.para[iaROC, :, 1]
    roc2  = ac2.para[iaROC, :, 1]
    fW1   = ac1.para[iafracW, :, 1]
    fW2   = ac2.para[iafracW, :, 1]
    M1    = ac1.para[iaMach, :, 1]
    M2    = ac2.para[iaMach, :, 1]
    gamV1 = ac1.para[iagamV, :, 1]
    gamV2 = ac2.para[iagamV, :, 1]

    sel_x = @lift Float64[$sel_pt]

    # Altitude
    ax1 = Axis(gl[1,1]; title="Altitude", ylabel="Alt [km]")
    lines!(ax1, pts, alt1; color=C1, linewidth=2.5, label=ac1.name)
    lines!(ax1, pts, alt2; color=C2, linewidth=2.5, label=ac2.name)
    vlines!(ax1, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)
    axislegend(ax1; position=:rb, labelsize=10)

    # Mach
    ax2 = Axis(gl[1,2]; title="Mach", ylabel="Mach")
    lines!(ax2, pts, M1; color=C1, linewidth=2.5)
    lines!(ax2, pts, M2; color=C2, linewidth=2.5)
    vlines!(ax2, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    # L/D
    ax3 = Axis(gl[2,1]; title="Lift-to-Drag", ylabel="L/D", xlabel="Mission Point")
    lines!(ax3, pts, LoD1; color=C1, linewidth=2.5)
    lines!(ax3, pts, LoD2; color=C2, linewidth=2.5)
    vlines!(ax3, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    # Rate of climb
    ax4 = Axis(gl[2,2]; title="Rate of Climb / Path Angle",
               ylabel="ROC [m/s]", xlabel="Mission Point")
    lines!(ax4, pts, roc1; color=C1, linewidth=2.5, label="ROC-1")
    lines!(ax4, pts, roc2; color=C2, linewidth=2.5, label="ROC-2")
    ax4r = Axis(gl[2,2]; ylabel="γ [rad]", yaxisposition=:right,
                backgroundcolor=:transparent)
    hidespines!(ax4r); hidexdecorations!(ax4r)
    lines!(ax4r, pts, gamV1; color=C1, linewidth=1.5, linestyle=:dot)
    lines!(ax4r, pts, gamV2; color=C2, linewidth=1.5, linestyle=:dot)
    hlines!(ax4, [0.0]; color=(RGBf(1,1,1),0.25), linewidth=1)
    vlines!(ax4, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    # Fuel fraction
    ax5 = Axis(gl[3,1]; title="Weight Fraction (W/WMTO)", ylabel="W/W₀", xlabel="Mission Point")
    lines!(ax5, pts, fW1; color=C1, linewidth=2.5)
    lines!(ax5, pts, fW2; color=C2, linewidth=2.5)
    vlines!(ax5, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    # CD×10⁴ along mission
    CD1 = ac1.para[iaCD, :, 1]; CD2 = ac2.para[iaCD, :, 1]
    ax6 = Axis(gl[3,2]; title="Total CD [×10⁻⁴]", ylabel="CD×10⁴", xlabel="Mission Point")
    lines!(ax6, pts, CD1.*1e4; color=C1, linewidth=2.5)
    lines!(ax6, pts, CD2.*1e4; color=C2, linewidth=2.5)
    vlines!(ax6, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    colgap!(gl, 8); rowgap!(gl, 8)
end

function build_engine_panel!(gl, ac1, ac2, sel_pt)
    # Station temperature arrays
    Tvars = [ieTt0, ieTt2, ieTt21, ieTt25, ieTt3, ieTt4, ieTt41, ieTt45, ieTt5, ieTt7]
    Tnames= ["Tt₀","Tt₂","Tt₂₁","Tt₂₅","Tt₃","Tt₄","Tt₄₁","Tt₄₅","Tt₅","Tt₇"]
    Pvars = [iept0, iept2, iept21, iept25, iept3, iept4, iept41, iept45, iept5, iept7]
    Pnames= ["pt₀","pt₂","pt₂₁","pt₂₅","pt₃","pt₄","pt₄₁","pt₄₅","pt₅","pt₇"]

    ns  = 1:length(Tvars)
    nsp = 1:length(Pvars)

    T1 = @lift [ac1.pare[v, $sel_pt, 1] for v in Tvars]
    T2 = @lift [ac2.pare[v, $sel_pt, 1] for v in Tvars]
    P1 = @lift [ac1.pare[v, $sel_pt, 1] / 1e3 for v in Pvars]   # kPa
    P2 = @lift [ac2.pare[v, $sel_pt, 1] / 1e3 for v in Pvars]

    ax_title = @lift "Engine Station Temperatures — $(ip_labels[$sel_pt])"

    # Temperature
    ax_T = Axis(gl[1, 1:2];
        title    = ax_title,
        xticks   = (collect(ns), Tnames),
        ylabel   = "Tt [K]",
        xticklabelrotation = 0.0)
    lines!(ax_T, collect(ns), T1; color=C1, linewidth=2.5, marker=:circle, markersize=9, label=ac1.name)
    lines!(ax_T, collect(ns), T2; color=C2, linewidth=2.5, marker=:circle, markersize=9, label=ac2.name)
    axislegend(ax_T; position=:lt, labelsize=10)

    # Pressure
    ax_P = Axis(gl[2, 1:2];
        title    = "Engine Station Total Pressures",
        xticks   = (collect(nsp), Pnames),
        ylabel   = "pt [kPa]")
    lines!(ax_P, collect(nsp), P1; color=C1, linewidth=2.5, marker=:circle, markersize=9)
    lines!(ax_P, collect(nsp), P2; color=C2, linewidth=2.5, marker=:circle, markersize=9)

    # Mission trends ──────────────────────────────────────────────────────────
    pts  = collect(1:iptotal)
    sel_x = @lift Float64[$sel_pt]

    Tt4_1 = ac1.pare[ieTt4,  :, 1]
    Tt4_2 = ac2.pare[ieTt4,  :, 1]
    BPR1  = ac1.pare[ieBPR,  :, 1]
    BPR2  = ac2.pare[ieBPR,  :, 1]
    OPR1  = ac1.pare[ieOPR,  :, 1]
    OPR2  = ac2.pare[ieOPR,  :, 1]
    TSFC1 = ac1.pare[ieTSFC, :, 1] .* 1e6
    TSFC2 = ac2.pare[ieTSFC, :, 1] .* 1e6
    FPR1  = ac1.pare[iepif,  :, 1]
    FPR2  = ac2.pare[iepif,  :, 1]
    mdot1 = ac1.pare[iemdotf, :, 1] .* 1e3   # g/s
    mdot2 = ac2.pare[iemdotf, :, 1] .* 1e3

    ax_tt4 = Axis(gl[3,1]; title="Tt4 Along Mission", ylabel="Tt4 [K]", xlabel="Mission Pt")
    lines!(ax_tt4, pts, Tt4_1; color=C1, linewidth=2.5, label=ac1.name)
    lines!(ax_tt4, pts, Tt4_2; color=C2, linewidth=2.5, label=ac2.name)
    vlines!(ax_tt4, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)
    axislegend(ax_tt4; position=:rt, labelsize=9)

    ax_tsfc = Axis(gl[3,2]; title="TSFC Along Mission", ylabel="TSFC [mg/N/s]", xlabel="Mission Pt")
    lines!(ax_tsfc, pts, TSFC1; color=C1, linewidth=2.5)
    lines!(ax_tsfc, pts, TSFC2; color=C2, linewidth=2.5)
    vlines!(ax_tsfc, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    ax_bpr = Axis(gl[4,1]; title="BPR + OPR/10 Along Mission", ylabel="BPR  /  OPR×0.1", xlabel="Mission Pt")
    lines!(ax_bpr, pts, BPR1;       color=C1, linewidth=2.5, label="BPR-1")
    lines!(ax_bpr, pts, BPR2;       color=C2, linewidth=2.5, label="BPR-2")
    lines!(ax_bpr, pts, OPR1.*0.1;  color=C1, linewidth=2.0, linestyle=:dash, label="OPR/10-1")
    lines!(ax_bpr, pts, OPR2.*0.1;  color=C2, linewidth=2.0, linestyle=:dash, label="OPR/10-2")
    vlines!(ax_bpr, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)
    axislegend(ax_bpr; position=:lb, labelsize=9, nbanks=2)

    ax_fpr = Axis(gl[4,2]; title="FPR & Fuel Flow Along Mission", ylabel="FPR / ṁ_fuel [g/s]", xlabel="Mission Pt")
    lines!(ax_fpr, pts, FPR1;   color=C1, linewidth=2.5, label="FPR-1")
    lines!(ax_fpr, pts, FPR2;   color=C2, linewidth=2.5, label="FPR-2")
    ax_fpr2 = Axis(gl[4,2]; ylabel="ṁ [g/s]", yaxisposition=:right, backgroundcolor=:transparent)
    hidespines!(ax_fpr2); hidexdecorations!(ax_fpr2)
    lines!(ax_fpr2, pts, mdot1; color=C1, linewidth=1.8, linestyle=:dot)
    lines!(ax_fpr2, pts, mdot2; color=C2, linewidth=1.8, linestyle=:dot)
    vlines!(ax_fpr, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    colgap!(gl, 8); rowgap!(gl, 8)
end

function build_drag_panel!(gl, ac1, ac2, sel_pt)
    drag_keys = ["CDi","CDfuse","CDwing","CDhtail","CDvtail","CDnace"]
    drag_inds = [iaCDi, iaCDfuse, iaCDwing, iaCDhtail, iaCDvtail, iaCDnace]
    drag_cols = [RGBf(0.27,0.51,0.71), RGBf(0.85,0.33,0.10),
                 RGBf(0.18,0.63,0.22), RGBf(0.95,0.60,0.07),
                 RGBf(0.60,0.20,0.70), RGBf(0.90,0.40,0.55)]
    S1 = ac1.wing.layout.S; S2 = ac2.wing.layout.S
    nx = length(drag_keys)

    # Reactive drag values at selected point (vector, one entry per component)
    cd1_obs = @lift [ac1.para[drag_inds[i], $sel_pt, 1] * S1 * 1e4 for i in 1:nx]
    cd2_obs = @lift [ac2.para[drag_inds[i], $sel_pt, 1] * S2 * 1e4 for i in 1:nx]

    ax_bar = Axis(gl[1, 1:2];
        title  = @lift("Drag Area (CD×S) at Mission Pt $($sel_pt): $(ip_labels[$sel_pt])"),
        xticks = (1:nx, drag_keys),
        ylabel = "CD×S [cm²]")

    # barplot! supports Observable y vectors; x positions are static
    barplot!(ax_bar, (1:nx) .- 0.22, cd1_obs;
        width=0.38, color=[(drag_cols[i],0.85) for i in 1:nx], label=ac1.name)
    barplot!(ax_bar, (1:nx) .+ 0.22, cd2_obs;
        width=0.38, color=[(drag_cols[i],0.55) for i in 1:nx],
        strokecolor=drag_cols, strokewidth=1, label=ac2.name)
    axislegend(ax_bar; position=:rt, labelsize=10)

    pts   = collect(1:iptotal)
    sel_x = @lift Float64[$sel_pt]

    # All mission points drag fractions (ac1)
    CD1_all = ac1.para[iaCD, :, 1]
    CD2_all = ac2.para[iaCD, :, 1]
    LoD1    = ac1.para[iaCL,:,1] ./ CD1_all
    LoD2    = ac2.para[iaCL,:,1] ./ CD2_all

    # Stacked drag fraction along mission (ac1)
    comp_fracs_1 = [ac1.para[idx, :, 1] ./ CD1_all for idx in drag_inds]
    comp_fracs_2 = [ac2.para[idx, :, 1] ./ CD2_all for idx in drag_inds]

    ax_stack1 = Axis(gl[2,1]; title="$(ac1.name) – Drag Fractions Along Mission",
                     xlabel="Mission Pt", ylabel="Fraction of CD")
    ax_stack2 = Axis(gl[2,2]; title="$(ac2.name) – Drag Fractions Along Mission",
                     xlabel="Mission Pt", ylabel="Fraction of CD")

    for (ax, fracs) in [(ax_stack1, comp_fracs_1), (ax_stack2, comp_fracs_2)]
        cumulative = zeros(length(pts))
        for (i, frac) in enumerate(fracs)
            band!(ax, pts, cumulative, cumulative .+ frac;
                color=(drag_cols[i], 0.70), label=drag_keys[i])
            cumulative .+= frac
        end
        vlines!(ax, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)
        axislegend(ax; position=:rb, labelsize=8, nbanks=2)
    end

    ax_LoD = Axis(gl[3,1]; title="L/D Along Mission", ylabel="L/D", xlabel="Mission Pt")
    lines!(ax_LoD, pts, LoD1; color=C1, linewidth=2.5, label=ac1.name)
    lines!(ax_LoD, pts, LoD2; color=C2, linewidth=2.5, label=ac2.name)
    vlines!(ax_LoD, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)
    axislegend(ax_LoD; position=:lb, labelsize=10)

    ax_CD = Axis(gl[3,2]; title="Total CD [×10⁻⁴] Along Mission",
                  ylabel="CD×10⁴", xlabel="Mission Pt")
    lines!(ax_CD, pts, CD1_all.*1e4; color=C1, linewidth=2.5)
    lines!(ax_CD, pts, CD2_all.*1e4; color=C2, linewidth=2.5)
    vlines!(ax_CD, sel_x; color=(YELLOW,0.8), linewidth=2, linestyle=:dash)

    colgap!(gl, 8); rowgap!(gl, 8)
end

function build_weights_panel!(gl, ac1, ac2)
    function get_weights(ac)
        pg = ac.parg; w = ac.wing; lg = ac.landing_gear
        Wbox  = w.inboard.webs.weight.W + w.inboard.caps.weight.W
        Wwing = Wbox * (1 + w.weight_frac_flap + w.weight_frac_slat +
                        w.weight_frac_ailerons + w.weight_frac_leading_trailing_edge +
                        w.weight_frac_ribs + w.weight_frac_spoilers +
                        w.weight_frac_attachments)
        Wadd = pg[igWMTO]*ac.fuselage.HPE_sys.W +
               lg.nose_gear.weight.W + lg.main_gear.weight.W
        lbls = ["Fuselage","Wing","H-Tail","V-Tail","Engine Sys","Fuel Tank","Misc"]
        vals = [ac.fuselage.weight, Wwing, ac.htail.weight, ac.vtail.weight,
                pg[igWtesys], pg[igWftank], Wadd] ./ 9.81e3
        return lbls, vals
    end

    wlabels, wv1 = get_weights(ac1)
    _,        wv2 = get_weights(ac2)
    nx = length(wlabels)
    w_cols = [RGBf(0.27,0.51,0.71), RGBf(0.18,0.63,0.22), RGBf(0.95,0.60,0.07),
              RGBf(0.60,0.20,0.70), RGBf(0.85,0.33,0.10), RGBf(0.90,0.40,0.55),
              RGBf(0.50,0.50,0.55)]

    ax_w = Axis(gl[1, 1:2];
        title  = "Structural Weight Breakdown [t]",
        xticks = (1:nx, wlabels),
        ylabel = "Weight [t]")

    for i in 1:nx
        barplot!(ax_w, [i-0.22], [wv1[i]]; width=0.38, color=(w_cols[i],0.85),
                 label=(i==1 ? ac1.name : nothing))
        barplot!(ax_w, [i+0.22], [wv2[i]]; width=0.38,
                 color=(w_cols[i],0.50), strokecolor=w_cols[i], strokewidth=1,
                 label=(i==1 ? ac2.name : nothing))
    end
    axislegend(ax_w; position=:rt, labelsize=10)

    # Pie charts ──────────────────────────────────────────────────────────────
    function draw_pie!(ax, ac, col)
        WMTO = ac.parg[igWMTO]; Wf = ac.parg[igWfuel]; Wp = ac.parg[igWpay]
        We = WMTO - Wf - Wp
        fracs  = [We, Wf, Wp] ./ WMTO
        labels = ["OEW", "Fuel", "Payload"]
        pcols  = [col, RGBf(0.90,0.40,0.10), RGBf(0.18,0.63,0.22)]
        hidedecorations!(ax); hidespines!(ax)
        cum = 0.0
        for (f, lab, pc) in zip(fracs, labels, pcols)
            θ1 = 2π*cum; θ2 = 2π*(cum+f)
            n  = max(Int(ceil(120*f)), 4)
            θs = LinRange(θ1, θ2, n)
            pts = vcat([Point2f(0,0)], [Point2f(cos(t), sin(t)) for t in θs])
            poly!(ax, pts; color=(pc,0.85), strokecolor=:white, strokewidth=1.5)
            θm = (θ1+θ2)/2
            text!(ax, 0.62*cos(θm), 0.62*sin(θm);
                text="$lab\n$(round(f*100,digits=1))%",
                fontsize=10, align=(:center,:center), color=:white)
            cum += f
        end
        text!(ax, 0.0, -1.35;
            text="MTOW = $(round(WMTO/9.81e3,digits=1)) t",
            fontsize=11, align=(:center,:center), color=TXT_COL)
        limits!(ax, -1.5, 1.5, -1.6, 1.4)
    end

    ax_pie1 = Axis(gl[2,1]; title=ac1.name, aspect=DataAspect())
    ax_pie2 = Axis(gl[2,2]; title=ac2.name, aspect=DataAspect())
    draw_pie!(ax_pie1, ac1, C1)
    draw_pie!(ax_pie2, ac2, C2)

    # Delta table ─────────────────────────────────────────────────────────────
    WMTO1=ac1.parg[igWMTO]; WMTO2=ac2.parg[igWMTO]
    metric_rows = [
        ("MTOW [t]",       WMTO1/9.81e3,    WMTO2/9.81e3),
        ("OEW [t]",       (WMTO1-ac1.parg[igWfuel]-ac1.parg[igWpay])/9.81e3,
                          (WMTO2-ac2.parg[igWfuel]-ac2.parg[igWpay])/9.81e3),
        ("Fuel [t]",       ac1.parg[igWfuel]/9.81e3,  ac2.parg[igWfuel]/9.81e3),
        ("Payload [t]",    ac1.parg[igWpay]/9.81e3,   ac2.parg[igWpay]/9.81e3),
        ("Wing [t]",       wv1[2], wv2[2]),
        ("Fuselage [t]",   wv1[1], wv2[1]),
    ]

    ax_tbl = Axis(gl[3, 1:2]; title="Key Weight Comparison (Δ = AC2 vs AC1)",
                   xticks=(1:length(metric_rows), [r[1] for r in metric_rows]),
                   ylabel="Weight [t]", xticklabelrotation=0.0)
    v1s = [r[2] for r in metric_rows]
    v2s = [r[3] for r in metric_rows]
    barplot!(ax_tbl, 1:length(v1s) .- 0.22, v1s; width=0.38, color=(C1,0.85), label=ac1.name)
    barplot!(ax_tbl, 1:length(v2s) .+ 0.22, v2s; width=0.38, color=(C2,0.85), label=ac2.name)
    axislegend(ax_tbl; position=:rt, labelsize=10)

    colgap!(gl, 8); rowgap!(gl, 8)
end

function build_wing_panel!(gl, ac1, ac2)
    function chord_dist(ac)
        w = ac.wing
        bo=w.layout.root_span; bs=w.layout.break_span; b=w.layout.span
        co=w.layout.root_chord; cs=co*w.inboard.λ; ct=co*w.outboard.λ
        return [0.0, bo/b, bs/b, 1.0], [co, co, cs, ct]
    end

    η1, c1 = chord_dist(ac1); η2, c2 = chord_dist(ac2)

    ax_chord = Axis(gl[1,1]; title="Wing Chord Distribution",
                    xlabel="η = 2y/b", ylabel="Chord [m]")
    lines!(ax_chord, η1, c1; color=C1, linewidth=2.5, marker=:circle, markersize=9, label=ac1.name)
    lines!(ax_chord, η2, c2; color=C2, linewidth=2.5, marker=:circle, markersize=9, label=ac2.name)
    axislegend(ax_chord; position=:tr, labelsize=10)

    # Fuselage cross-section
    ax_fuse = Axis(gl[1,2]; title="Fuselage Cross-Section",
                   xlabel="y [m]", ylabel="z [m]", aspect=DataAspect())
    R1=ac1.fuselage.layout.radius; R2=ac2.fuselage.layout.radius
    θ = LinRange(0, 2π, 200)
    lines!(ax_fuse, R1.*cos.(θ), R1.*sin.(θ); color=C1, linewidth=2.5,
           label="$(ac1.name) R=$(round(R1,digits=3))m")
    lines!(ax_fuse, R2.*cos.(θ), R2.*sin.(θ); color=C2, linewidth=2.5,
           linestyle=:dash, label="$(ac2.name) R=$(round(R2,digits=3))m")
    axislegend(ax_fuse; position=:lb, labelsize=9)

    # Geometry comparison bar
    Lf1=ac1.fuselage.layout.x_end-ac1.fuselage.layout.x_nose
    Lf2=ac2.fuselage.layout.x_end-ac2.fuselage.layout.x_nose
    glabels=["Span [m]","Fuse L [m]","S [m²]","Root c [m]","AR","Sweep [°]"]
    gv1=[ac1.wing.layout.span, Lf1, ac1.wing.layout.S,
         ac1.wing.layout.root_chord, ac1.wing.layout.AR, ac1.wing.layout.sweep]
    gv2=[ac2.wing.layout.span, Lf2, ac2.wing.layout.S,
         ac2.wing.layout.root_chord, ac2.wing.layout.AR, ac2.wing.layout.sweep]

    ax_geom = Axis(gl[2, 1:2];
        title="Key Geometry Parameters",
        xticks=(1:length(glabels), glabels),
        ylabel="Value [m / m² / –]")
    barplot!(ax_geom, 1:length(gv1) .- 0.22, gv1; width=0.38, color=(C1,0.85), label=ac1.name)
    barplot!(ax_geom, 1:length(gv2) .+ 0.22, gv2; width=0.38, color=(C2,0.85), label=ac2.name)
    axislegend(ax_geom; position=:tr, labelsize=10)

    # Stability bar
    function SM_data(ac)
        xNP=ac.parg[igxNP]; xCGf=ac.parg[igxCGfwd]; xCGa=ac.parg[igxCGaft]
        mac=ac.wing.mean_aero_chord
        [(xNP-xCGf)/mac*100, (xNP-xCGa)/mac*100]
    end
    sm1=SM_data(ac1); sm2=SM_data(ac2)

    ax_sm = Axis(gl[3, 1:2];
        title="Static Margin [%MAC]  (fwd / aft CG)",
        xticks=(1:2, ["SM fwd CG","SM aft CG"]),
        ylabel="SM [%MAC]")
    barplot!(ax_sm, [1,2] .- 0.22, sm1; width=0.38, color=(C1,0.85), label=ac1.name)
    barplot!(ax_sm, [1,2] .+ 0.22, sm2; width=0.38, color=(C2,0.85), label=ac2.name)
    axislegend(ax_sm; position=:rt, labelsize=10)

    colgap!(gl, 8); rowgap!(gl, 8)
end

# ─────────────────────────────────────────────────────────────────────────────
# Top-level overview text panel
# ─────────────────────────────────────────────────────────────────────────────

function build_overview_panel!(gl, ac1, ac2)
    """
    A scrollable or fixed text comparison table using a simple Axis with
    text! calls (since Makie does not have a native DataTable widget).
    """
    ax = Axis(gl[1,1:2]; title="Key Performance & Geometry Comparison")
    hidedecorations!(ax); hidespines!(ax)
    limits!(ax, 0, 1, 0, 1)

    function fmt(v)
        abs(v) < 1e-10 && return "—"
        @sprintf("%.4g", v)
    end
    function Δ(v1, v2)
        v1 == 0 && return "—"
        d = (v2-v1)/abs(v1)*100
        d > 0 && return @sprintf("+%.1f%%", d)
        return @sprintf("%.1f%%", d)
    end

    para1 = view(ac1.para, :,:,1); pare1 = view(ac1.pare, :,:,1)
    para2 = view(ac2.para, :,:,1); pare2 = view(ac2.pare, :,:,1)

    rows = [
        ("",                   "Parameter",                      ac1.name,                              ac2.name,                              "Δ (2 vs 1)"),
        ("Performance",        "PFEI [kJ/kg/km]",               fmt(ac1.parm[imPFEI,1]),               fmt(ac2.parm[imPFEI,1]),               Δ(ac1.parm[imPFEI,1],   ac2.parm[imPFEI,1])),
        ("",                   "Design Range [km]",             fmt(ac1.parm[imRange,1]/1e3),          fmt(ac2.parm[imRange,1]/1e3),          Δ(ac1.parm[imRange,1],  ac2.parm[imRange,1])),
        ("",                   "Cruise L/D",                    fmt(para1[iaCL,ipcruise1]/para1[iaCD,ipcruise1]), fmt(para2[iaCL,ipcruise1]/para2[iaCD,ipcruise1]), Δ(para1[iaCL,ipcruise1]/para1[iaCD,ipcruise1], para2[iaCL,ipcruise1]/para2[iaCD,ipcruise1])),
        ("",                   "Cruise Mach",                   fmt(para1[iaMach,ipcruise1]),          fmt(para2[iaMach,ipcruise1]),          Δ(para1[iaMach,ipcruise1], para2[iaMach,ipcruise1])),
        ("Weights",            "MTOW [t]",                      fmt(ac1.parg[igWMTO]/9.81e3),          fmt(ac2.parg[igWMTO]/9.81e3),          Δ(ac1.parg[igWMTO],     ac2.parg[igWMTO])),
        ("",                   "OEW [t]",                       fmt((ac1.parg[igWMTO]-ac1.parg[igWfuel]-ac1.parg[igWpay])/9.81e3), fmt((ac2.parg[igWMTO]-ac2.parg[igWfuel]-ac2.parg[igWpay])/9.81e3), Δ(ac1.parg[igWMTO]-ac1.parg[igWfuel]-ac1.parg[igWpay], ac2.parg[igWMTO]-ac2.parg[igWfuel]-ac2.parg[igWpay])),
        ("",                   "Fuel [t]",                      fmt(ac1.parg[igWfuel]/9.81e3),         fmt(ac2.parg[igWfuel]/9.81e3),         Δ(ac1.parg[igWfuel],    ac2.parg[igWfuel])),
        ("Wing",               "Span [m]",                      fmt(ac1.wing.layout.span),             fmt(ac2.wing.layout.span),             Δ(ac1.wing.layout.span, ac2.wing.layout.span)),
        ("",                   "AR",                            fmt(ac1.wing.layout.AR),               fmt(ac2.wing.layout.AR),               Δ(ac1.wing.layout.AR,   ac2.wing.layout.AR)),
        ("",                   "Sweep [°]",                     fmt(ac1.wing.layout.sweep),            fmt(ac2.wing.layout.sweep),            Δ(ac1.wing.layout.sweep,ac2.wing.layout.sweep)),
        ("",                   "S [m²]",                        fmt(ac1.wing.layout.S),                fmt(ac2.wing.layout.S),                Δ(ac1.wing.layout.S,    ac2.wing.layout.S)),
        ("Engine (cruise)",    "Fan dia [m]",                   fmt(ac1.parg[igdfan]),                 fmt(ac2.parg[igdfan]),                 Δ(ac1.parg[igdfan],     ac2.parg[igdfan])),
        ("",                   "BPR",                           fmt(pare1[ieBPR,  ipcruise1]),         fmt(pare2[ieBPR,  ipcruise1]),         Δ(pare1[ieBPR,ipcruise1],  pare2[ieBPR,ipcruise1])),
        ("",                   "OPR",                           fmt(pare1[ieOPR,  ipcruise1]),         fmt(pare2[ieOPR,  ipcruise1]),         Δ(pare1[ieOPR,ipcruise1],  pare2[ieOPR,ipcruise1])),
        ("",                   "Tt4 [K]",                       fmt(pare1[ieTt4,  ipcruise1]),         fmt(pare2[ieTt4,  ipcruise1]),         Δ(pare1[ieTt4,ipcruise1],  pare2[ieTt4,ipcruise1])),
        ("",                   "TSFC [mg/N/s]",                 fmt(pare1[ieTSFC, ipcruise1]*1e6),     fmt(pare2[ieTSFC, ipcruise1]*1e6),     Δ(pare1[ieTSFC,ipcruise1], pare2[ieTSFC,ipcruise1])),
    ]

    ncols   = 4   # parameter | ac1 | ac2 | Δ
    col_x   = [0.01, 0.30, 0.55, 0.78]
    nrow    = length(rows)
    row_h   = 1.0 / (nrow+1)

    for (ri, row) in enumerate(rows)
        y = 1.0 - ri * row_h
        # group label (col 0 if non-empty)
        if !isempty(row[1])
            text!(ax, 0.0, y + row_h*0.2; text=row[1],
                  fontsize=10, color=YELLOW, font=:bold, align=(:left,:center))
        end
        for (ci, (x, val)) in enumerate(zip(col_x, row[2:end]))
            is_hdr = ri == 1
            col  = is_hdr ? TXT_COL : (ci==4 && val!="—" && startswith(val,"+") ? RGBf(0.3,0.85,0.3) :
                                        ci==4 && val!="—" && !startswith(val,"+") ? RGBf(0.9,0.35,0.35) : TXT_COL)
            text!(ax, x, y; text=val,
                  fontsize= is_hdr ? 11 : 10,
                  font    = is_hdr ? :bold : :regular,
                  color   = col,
                  align   = (:left,:center))
        end
        if ri == 1 || mod(ri, 2) == 0
            poly!(ax, [Point2f(0, y-row_h*0.05), Point2f(1,y-row_h*0.05),
                       Point2f(1, y+row_h*0.95), Point2f(0,y+row_h*0.95)];
                  color=(RGBf(1,1,1), ri==1 ? 0.06 : 0.03), strokewidth=0)
        end
    end

    colgap!(gl, 8)
end

# ─────────────────────────────────────────────────────────────────────────────
# Main layout builder
# ─────────────────────────────────────────────────────────────────────────────

"""
    launch_dashboard(ac1, ac2)

Open an interactive GLMakie comparison window for two sized TASOPT aircraft.
Returns the Figure object (window stays open until closed by the user).
"""
function launch_dashboard(ac1::TASOPT.aircraft, ac2::TASOPT.aircraft)

    @assert ac1.is_sized[1] "Aircraft 1 must be sized before comparison"
    @assert ac2.is_sized[1] "Aircraft 2 must be sized before comparison"

    # ── Reactive state ────────────────────────────────────────────────────────
    sel_pt     = Observable(ipcruise1)   # active mission point
    active_tab = Observable(:overview)  # current content tab

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = Figure(size=(1540, 980), backgroundcolor=DARK_BG)

    # Row 1 – title
    Label(fig[1, 1:2];
        text="✈  TASOPT Interactive Comparison Dashboard",
        fontsize=17, color=TXT_COL, font=:bold, tellwidth=false)
    Label(fig[1, 3];
        text="$(ac1.name)  vs  $(ac2.name)",
        fontsize=12, color=RGBf(0.7,0.7,0.7), tellwidth=false)

    # Row 2 – aircraft outlines (columns 1-2) + usage hint (column 3)
    ax_ac1 = Axis(fig[2, 1];
        title         = ac1.name,
        titlecolor    = C1,
        titlesize     = 13,
        xlabel        = "x [m]", ylabel="y [m]",
        aspect        = DataAspect(),
        backgroundcolor = AX_BG)
    ax_ac2 = Axis(fig[2, 2];
        title         = ac2.name,
        titlecolor    = C2,
        titlesize     = 13,
        xlabel        = "x [m]", ylabel="y [m]",
        aspect        = DataAspect(),
        backgroundcolor = AX_BG)

    comps1 = draw_aircraft!(ax_ac1, ac1, C1)
    comps2 = draw_aircraft!(ax_ac2, ac2, C2)

    # Usage hint panel
    hint_ax = Axis(fig[2, 3]; backgroundcolor=MID_BG)
    hidedecorations!(hint_ax); hidespines!(hint_ax)
    limits!(hint_ax, 0, 1, 0, 1)
    hint_lines = [
        "Hover → component tooltip",
        "",
        "Click wing  → Drag tab",
        "Click fuse  → Weights tab",
        "Click engine → Engine tab",
        "Click h-tail → Wing tab",
        "",
        "Slider → select mission pt",
        "Tabs   → switch analysis",
        "",
        "Yellow ▼ = Neutral Point",
        "White ─┼─ = CG range",
    ]
    for (i, ln) in enumerate(hint_lines)
        text!(hint_ax, 0.07, 0.93 - (i-1)*0.072;
            text=ln, fontsize=11, color=TXT_COL, align=(:left,:top))
    end

    # ── Enable DataInspector (hover tooltips on all poly! objects) ─────────────
    DataInspector(fig;
        indicator_color    = YELLOW,
        indicator_linewidth = 2.5,
        fontsize           = 11,
        backgroundcolor    = RGBAf(0.1,0.1,0.13,0.92),
        textcolor          = TXT_COL)

    # ── Click on outline components → switch tab ─────────────────────────────
    for (ax, comps) in [(ax_ac1, comps1), (ax_ac2, comps2)]
        on(events(ax.scene).mousebutton) do event
            event.button  == Mouse.left   || return
            event.action  == Mouse.press  || return
            mpos  = events(ax.scene).mouseposition[]
            hit, _ = pick(ax.scene, mpos)
            hit === nothing && return
            if hit === comps.wing
                active_tab[] = :drag
            elseif hit === comps.fuselage
                active_tab[] = :weights
            elseif hit === comps.htail
                active_tab[] = :wing
            elseif hit ∈ comps.engines
                active_tab[] = :engine
            end
        end
    end

    # ── Row 3 – mission-point slider ─────────────────────────────────────────
    sl_grid = fig[3, 1:3] = GridLayout()
    Label(sl_grid[1,1]; text="Mission Point:", color=TXT_COL, fontsize=12, tellwidth=false)
    sl = Slider(sl_grid[1,2]; range=1:iptotal, startvalue=ipcruise1,
                color_active=RGBf(0.27,0.51,0.71),
                color_inactive=RGBf(0.30,0.30,0.35))
    pt_lbl = @lift "$(ip_labels[$sel_pt])  (pt $($sel_pt))"
    Label(sl_grid[1,3]; text=pt_lbl, color=YELLOW, fontsize=12, tellwidth=false)
    on(sl.value) do v; sel_pt[] = v; end
    colsize!(sl_grid, 2, Relative(0.65))

    # ── Row 4 – tab buttons ───────────────────────────────────────────────────
    tab_defs = [
        (:overview, "📊 Overview"),
        (:mission,  "✈ Mission Profile"),
        (:engine,   "⚙ Engine Analysis"),
        (:drag,     "💨 Drag Breakdown"),
        (:weights,  "⚖ Weights"),
        (:wing,     "📐 Wing & Geometry"),
    ]
    btn_grid = fig[4, 1:3] = GridLayout()
    for (i, (sym, lbl)) in enumerate(tab_defs)
        active_col = @lift $active_tab == sym ? RGBf(0.27,0.51,0.71) : RGBf(0.28,0.28,0.33)
        btn = Button(btn_grid[1,i];
            label          = lbl,
            buttoncolor    = active_col,
            labelcolor     = TXT_COL,
            fontsize       = 11,
            padding        = (8,8,5,5))
        on(btn.clicks) do _; active_tab[] = sym; end
    end

    # ── Row 5 – dynamic content area ─────────────────────────────────────────
    content_gl = fig[5, 1:3] = GridLayout()

    function rebuild!(tab)
        # Delete all Axis objects attached to the content GridLayout, then empty it.
        # Axes are disposed properly; GLMakie GC handles the rest.
        for item in copy(content_gl.content)
            ob = item.content
            try
                if ob isa Axis
                    delete!(ob)
                end
            catch
            end
        end
        empty!(content_gl)

        if     tab == :overview;  build_overview_panel!(content_gl, ac1, ac2)
        elseif tab == :mission;   build_mission_panel!( content_gl, ac1, ac2, sel_pt)
        elseif tab == :engine;    build_engine_panel!(  content_gl, ac1, ac2, sel_pt)
        elseif tab == :drag;      build_drag_panel!(    content_gl, ac1, ac2, sel_pt)
        elseif tab == :weights;   build_weights_panel!( content_gl, ac1, ac2)
        elseif tab == :wing;      build_wing_panel!(    content_gl, ac1, ac2)
        end
    end

    on(active_tab) do tab; rebuild!(tab); end
    rebuild!(active_tab[])   # initial build

    # ── Row/col sizing ────────────────────────────────────────────────────────
    rowsize!(fig.layout, 1, Auto(0.04))   # title
    rowsize!(fig.layout, 2, Auto(0.30))   # aircraft outlines
    rowsize!(fig.layout, 3, Auto(0.04))   # slider
    rowsize!(fig.layout, 4, Auto(0.04))   # tab buttons
    rowsize!(fig.layout, 5, Auto(0.58))   # content

    colsize!(fig.layout, 1, Relative(0.35))
    colsize!(fig.layout, 2, Relative(0.35))
    colsize!(fig.layout, 3, Relative(0.30))

    display(fig)
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

function main()
    if length(ARGS) < 2
        println("""
        TASOPT Interactive Comparison Dashboard

        Usage:
            julia --project=.. interactive_dashboard.jl <ac1.jld2> <ac2.jld2>

        Save a sized aircraft with:
            TASOPT.quicksave_aircraft(ac, "myfile.jld2")

        From the REPL:
            include("interactive_dashboard.jl")
            launch_dashboard(ac1, ac2)
        """)
        return
    end

    ac1_path, ac2_path = ARGS[1], ARGS[2]
    println("Loading $ac1_path …"); ac1 = TASOPT.quickload_aircraft(ac1_path)
    println("Loading $ac2_path …"); ac2 = TASOPT.quickload_aircraft(ac2_path)

    println("Launching dashboard …")
    fig = launch_dashboard(ac1, ac2)
    wait(display(fig))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
