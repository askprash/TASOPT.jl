# compare_aircraft_script.jl
#
# Standalone (non-Pluto) script to generate a multi-panel comparison
# PDF/PNG between two sized TASOPT aircraft saved as .jld2 files.
#
# Usage:
#   julia --project=.. compare_aircraft_script.jl path/to/ac1.jld2 path/to/ac2.jld2
#
# Or from the REPL:
#   include("compare_aircraft_script.jl")
#   compare_and_save(ac1, ac2; output="comparison.pdf")

using TASOPT
using Plots, StatsPlots, Printf
include(TASOPT.__TASOPTindices__)

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

function extract_weights(ac)
    parg = ac.parg
    fuse = ac.fuselage; wing = ac.wing; htail = ac.htail; vtail = ac.vtail
    lg   = ac.landing_gear

    Wbox  = wing.inboard.webs.weight.W + wing.inboard.caps.weight.W
    Wwing = Wbox * (1 + wing.weight_frac_flap + wing.weight_frac_slat +
                    wing.weight_frac_ailerons + wing.weight_frac_leading_trailing_edge +
                    wing.weight_frac_ribs + wing.weight_frac_spoilers +
                    wing.weight_frac_attachments)
    Wadd = parg[igWMTO]*fuse.HPE_sys.W + lg.nose_gear.weight.W + lg.main_gear.weight.W

    labels = ["Fuselage", "Wing", "H-Tail", "V-Tail", "Engine Sys", "Fuel Tank", "Added"]
    vals   = [fuse.weight, Wwing, htail.weight, vtail.weight,
              parg[igWtesys], parg[igWftank], Wadd]
    return labels, vals
end

function wing_chord_dist(ac)
    w  = ac.wing
    bo = w.layout.root_span; bs = w.layout.break_span; b = w.layout.span
    co = w.layout.root_chord
    cs = co * w.inboard.λ;   ct = co * w.outboard.λ
    η  = [0, bo/b, bs/b, 1.0]
    c  = [co, co, cs, ct]
    return η, c
end

# -----------------------------------------------------------------------
# Main comparison function
# -----------------------------------------------------------------------

"""
    compare_and_save(ac1, ac2; output="comparison.pdf")

Generate a multi-panel comparison figure for two sized aircraft and save it.
"""
function compare_and_save(ac1::TASOPT.aircraft, ac2::TASOPT.aircraft;
                          output::String = "comparison.pdf")

    c1 = :steelblue
    c2 = :tomato

    # ── 1. Stick figures ──────────────────────────────────────────────
    p_sf1 = TASOPT.stickfig(ac1; annotate_text=true,  annotate_length=true,
                                  annotate_group=false, label_fs=9)
    title!(p_sf1, ac1.name)

    p_sf2 = TASOPT.stickfig(ac2; annotate_text=true,  annotate_length=true,
                                  annotate_group=false, label_fs=9)
    title!(p_sf2, ac2.name)

    # ── 2. Weight breakdown ────────────────────────────────────────────
    wlabels, wv1 = extract_weights(ac1)
    _,        wv2 = extract_weights(ac2)

    p_w = groupedbar(wlabels, [wv1 wv2] ./ 9.81e3;
        bar_position=:dodge, bar_width=0.7,
        label=[ac1.name ac2.name], color=[c1 c2], alpha=0.85,
        ylabel="Weight [t]", title="Weight Breakdown",
        legend=:topright, grid=:y, gridcolor=:lightgray)

    # ── 3. Drag breakdown ──────────────────────────────────────────────
    drag_keys = ["CDi", "CDfuse", "CDwing", "CDhtail", "CDvtail", "CDnace"]
    drag_inds = [iaCDi, iaCDfuse, iaCDwing, iaCDhtail, iaCDvtail, iaCDnace]

    cd1 = [ac1.para[i, ipcruise1, 1] for i in drag_inds] .* ac1.wing.layout.S .* 1e4
    cd2 = [ac2.para[i, ipcruise1, 1] for i in drag_inds] .* ac2.wing.layout.S .* 1e4

    p_cd = groupedbar(drag_keys, [cd1 cd2];
        bar_position=:dodge, bar_width=0.7,
        label=[ac1.name ac2.name], color=[c1 c2], alpha=0.85,
        ylabel="CD×S [cm²]", title="Cruise Drag Areas",
        legend=:topright, grid=:y, gridcolor=:lightgray)

    # ── 4. Mission profile ─────────────────────────────────────────────
    npts = size(ac1.para, 2)
    pts  = 1:npts

    alt1  = ac1.para[iaalt, :, 1] ./ 1e3
    alt2  = ac2.para[iaalt, :, 1] ./ 1e3
    LoD1  = ac1.para[iaCL, :, 1] ./ ac1.para[iaCD, :, 1]
    LoD2  = ac2.para[iaCL, :, 1] ./ ac2.para[iaCD, :, 1]

    p_alt = plot(pts, alt1; label=ac1.name, color=c1, lw=2,
                 ylabel="Altitude [km]", title="Altitude Profile", legend=:bottomright)
    plot!(p_alt, pts, alt2; label=ac2.name, color=c2, lw=2)

    p_LoD = plot(pts, LoD1; label=ac1.name, color=c1, lw=2,
                 ylabel="L/D", xlabel="Flight Point", title="L/D", legend=:bottomright)
    plot!(p_LoD, pts, LoD2; label=ac2.name, color=c2, lw=2)

    # ── 5. Wing chord distribution ─────────────────────────────────────
    η1, cw1 = wing_chord_dist(ac1)
    η2, cw2 = wing_chord_dist(ac2)

    p_chord = plot(η1, cw1; label=ac1.name, color=c1, lw=2.5, marker=:circle, ms=5,
                   xlabel="η", ylabel="Chord [m]", title="Wing Chord Distribution",
                   legend=:topright)
    plot!(p_chord, η2, cw2; label=ac2.name, color=c2, lw=2.5, marker=:circle, ms=5)

    # ── 6. Key scalar comparison bar ──────────────────────────────────
    scalar_labels = ["PFEI\n[kJ/kg/km]", "L/D\ncruise", "TSFC×100\n[mg/N/s]",
                     "BPR",              "OPR/10"]
    sv1 = [ac1.parm[imPFEI,1],
           ac1.para[iaCL,ipcruise1,1]/ac1.para[iaCD,ipcruise1,1],
           ac1.pare[ieTSFC,ipcruise1,1]*1e8,
           ac1.pare[ieBPR, ipcruise1,1],
           ac1.pare[ieOPR, ipcruise1,1]/10]
    sv2 = [ac2.parm[imPFEI,1],
           ac2.para[iaCL,ipcruise1,1]/ac2.para[iaCD,ipcruise1,1],
           ac2.pare[ieTSFC,ipcruise1,1]*1e8,
           ac2.pare[ieBPR, ipcruise1,1],
           ac2.pare[ieOPR, ipcruise1,1]/10]

    p_scalar = groupedbar(scalar_labels, [sv1 sv2];
        bar_position=:dodge, bar_width=0.7,
        label=[ac1.name ac2.name], color=[c1 c2], alpha=0.85,
        ylabel="Value (scaled)", title="Key Performance Indicators",
        legend=:topright, grid=:y, gridcolor=:lightgray)

    # ── Assemble composite figure ──────────────────────────────────────
    fig = plot(
        p_sf1,   p_sf2,
        p_w,     p_cd,
        p_alt,   p_LoD,
        p_chord, p_scalar;
        layout=(4, 2),
        size=(1200, 1600),
        dpi=150,
        leftmargin=6Plots.mm,
        bottommargin=5Plots.mm,
        titlelocation=:center
    )

    savefig(fig, output)
    println("Comparison saved to: $output")
    return fig
end

# -----------------------------------------------------------------------
# CLI entry point
# -----------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: julia --project=.. compare_aircraft_script.jl <ac1.jld2> <ac2.jld2> [output.pdf]")
        exit(1)
    end
    ac1_path = ARGS[1]
    ac2_path = ARGS[2]
    out_path = length(ARGS) >= 3 ? ARGS[3] : "comparison.pdf"

    println("Loading $(ac1_path) ...")
    ac1 = TASOPT.quickload_aircraft(ac1_path)
    println("Loading $(ac2_path) ...")
    ac2 = TASOPT.quickload_aircraft(ac2_path)

    compare_and_save(ac1, ac2; output=out_path)
end
