export export_dashboard_data, export_dashboard_bundle, start_dashboard

const DASHBOARD_SCHEMA_VERSION = "1.5.0"
const DASHBOARD_PLOTLY_SRC = "https://cdn.plot.ly/plotly-2.35.0.min.js"

const DASHBOARD_MISSION_LABELS = [
    "Static", "Rotate", "Takeoff", "Cutback",
    "Climb 1", "Climb 2", "Climb 3", "Climb 4", "Climb 5",
    "Cruise 1", "Cruise 2",
    "Descent 1", "Descent 2", "Descent 3", "Descent 4", "Descent 5",
    "Test",
]

const DASHBOARD_MISSION_SHORT = [
    "ST", "RO", "TO", "CB", "CL1", "CL2", "CL3", "CL4", "CL5",
    "CR1", "CR2", "DS1", "DS2", "DS3", "DS4", "DS5", "TST",
]

const DASHBOARD_DRAG_KEYS = ["CDi", "CDfuse", "CDwing", "CDhtail", "CDvtail", "CDnace"]
const DASHBOARD_DRAG_INDICES = [iaCDi, iaCDfuse, iaCDwing, iaCDhtail, iaCDvtail, iaCDnace]

const DASHBOARD_AEROPERF_FIELDS = [
    (key = :LDs, label = "L/D"),
    (key = :CDs, label = "Total CD"),
    (key = :CLhs, label = "Tail CLh"),
    (key = :CDis, label = "Induced drag CDi"),
    (key = :CDwings, label = "Wing drag CD"),
    (key = :CDfuses, label = "Fuselage drag CD"),
    (key = :CDhtails, label = "H-tail drag CD"),
    (key = :CDvtails, label = "V-tail drag CD"),
    (key = :CDothers, label = "Other drag CD"),
    (key = :clpos, label = "Section cl root"),
    (key = :clpss, label = "Section cl break"),
    (key = :clpts, label = "Section cl tip"),
    (key = :cdfss, label = "Section cd friction"),
    (key = :cdpss, label = "Section cd pressure"),
    (key = :cdwss, label = "Section cd wave"),
    (key = :cdss, label = "Section cd total"),
]

const DASHBOARD_T_STATION_LABELS = ["Tt0", "Tt2", "Tt21", "Tt25", "Tt3", "Tt4", "Tt41", "Tt45", "Tt5", "Tt7"]
const DASHBOARD_T_STATION_INDICES = [ieTt0, ieTt2, ieTt21, ieTt25, ieTt3, ieTt4, ieTt41, ieTt45, ieTt5, ieTt7]

const DASHBOARD_P_STATION_LABELS = ["Pt0", "Pt2", "Pt21", "Pt25", "Pt3", "Pt4", "Pt41", "Pt45", "Pt5", "Pt7"]
const DASHBOARD_P_STATION_INDICES = [iept0, iept2, iept21, iept25, iept3, iept4, iept41, iept45, iept5, iept7]

safe_float(x) = (x isa Real && isfinite(x)) ? Float64(x) : 0.0
safe_array(x) = [safe_float(v) for v in vec(x)]
tonnes_from_force(x) = safe_float(x / 9.81e3)
kg_per_ns_to_mg_per_ns(x) = safe_float(x * 1e6)
safe_ratio(num, den) = abs(den) < eps(Float64) ? 0.0 : safe_float(num / den)
safe_nested(x::AbstractMatrix) = [[begin
    value = x[i, j]
    (value isa Real && isfinite(value)) ? Float64(value) : missing
end for j in axes(x, 2)] for i in axes(x, 1)]

first_column_vector(x) = ndims(x) == 1 ? vec(x) : vec(view(x, :, 1))
mission_matrix(x) = ndims(x) == 2 ? x : view(x, :, :, 1)

function raw_index_entries(prefix::AbstractString, total::Integer)
    total_symbol = Symbol(prefix * "total")
    pairs = Tuple{Int, String}[]
    seen = Set{Int}()

    for sym in names(@__MODULE__; all = true, imported = false)
        sym == total_symbol && continue
        name = String(sym)
        startswith(name, prefix) || continue
        isdefined(@__MODULE__, sym) || continue
        value = getfield(@__MODULE__, sym)
        value isa Integer || continue
        index = Int(value)
        1 <= index <= total || continue
        index in seen && continue
        push!(pairs, (index, name))
        push!(seen, index)
    end

    sort!(pairs; by = first)
    return pairs
end

function raw_scalar_entries(values, prefix::AbstractString, source::AbstractString)
    entries = Any[]
    data = safe_array(values)
    for (index, key) in raw_index_entries(prefix, length(data))
        push!(entries, (
            key = key,
            source = String(source),
            index = index,
            value = data[index],
        ))
    end
    return entries
end

function raw_series_entries(values, prefix::AbstractString, source::AbstractString)
    entries = Any[]
    matrix = mission_matrix(values)
    row_count = size(matrix, 1)
    for (index, key) in raw_index_entries(prefix, row_count)
        push!(entries, (
            key = key,
            source = String(source),
            index = index,
            values = safe_array(view(matrix, index, :)),
        ))
    end
    return entries
end

function raw_payload(ac::aircraft)
    return (
        note = "Raw labeled TASOPT variables for advanced dashboard plots. ig*/im* are scalar values; ia*/ie* are mission-point series.",
        scalar_values = vcat(
            raw_scalar_entries(ac.parg, "ig", "parg"),
            raw_scalar_entries(first_column_vector(ac.parm), "im", "parm"),
        ),
        mission_series = vcat(
            raw_series_entries(ac.para, "ia", "para"),
            raw_series_entries(ac.pare, "ie", "pare"),
        ),
    )
end

function json_escape(s::AbstractString)
    escaped = escape_string(String(s))
    escaped = replace(escaped, "<" => "\\u003c", ">" => "\\u003e", "&" => "\\u0026")
    return escaped
end

json_indent(level::Integer) = repeat(" ", level)
json_scalar_like(x) = x isa Union{String, Bool, Nothing, Missing, Real, Symbol}

_to_json(x) = _to_json(x, 0)
_to_json(x::String, ::Integer) = "\"" * json_escape(x) * "\""
_to_json(x::Bool, ::Integer) = x ? "true" : "false"
_to_json(x::Nothing, ::Integer) = "null"
_to_json(x::Missing, ::Integer) = "null"
_to_json(x::Real, ::Integer) = isnan(x) || isinf(x) ? "null" : string(x)
_to_json(x::Symbol, indent::Integer) = _to_json(String(x), indent)

function _to_json(items::AbstractVector, indent::Integer)
    isempty(items) && return "[]"

    if all(json_scalar_like, items)
        return "[" * join((_to_json(item, indent) for item in items), ", ") * "]"
    end

    inner_indent = indent + 2
    body = join((json_indent(inner_indent) * _to_json(item, inner_indent) for item in items), ",\n")
    return "[\n" * body * "\n" * json_indent(indent) * "]"
end

_to_json(items::Tuple, indent::Integer) = _to_json(collect(items), indent)

function _to_json(d::AbstractDict, indent::Integer)
    isempty(d) && return "{}"
    inner_indent = indent + 2
    body = join((
        json_indent(inner_indent) * _to_json(string(k), inner_indent) * ": " * _to_json(v, inner_indent)
        for (k, v) in d
    ), ",\n")
    return "{\n" * body * "\n" * json_indent(indent) * "}"
end

function _to_json(nt::NamedTuple, indent::Integer)
    isempty(nt) && return "{}"
    inner_indent = indent + 2
    body = join((
        json_indent(inner_indent) * _to_json(string(k), inner_indent) * ": " * _to_json(v, inner_indent)
        for (k, v) in pairs(nt)
    ), ",\n")
    return "{\n" * body * "\n" * json_indent(indent) * "}"
end

function dashboard_dist_dir()
    return joinpath(dirname(__TASOPTroot__), "dashboard", "dist")
end

function ensure_dashboard_dist()
    dist_dir = dashboard_dist_dir()
    assets_dir = joinpath(dist_dir, "assets")
    isfile(joinpath(dist_dir, "index.html")) || throw(ArgumentError(
        "Dashboard build output is missing at $(joinpath(dist_dir, "index.html")). Run `cd dashboard && npm run build` first."))
    isdir(assets_dir) || throw(ArgumentError(
        "Dashboard asset directory is missing at $assets_dir. Run `cd dashboard && npm run build` first."))
    return dist_dir, assets_dir
end

function copy_dashboard_dist!(target_dir::AbstractString)
    dist_dir, _ = ensure_dashboard_dist()
    mkpath(target_dir)
    for (root, _, files) in walkdir(dist_dir)
        rel_root = relpath(root, dist_dir)
        target_root = rel_root == "." ? target_dir : joinpath(target_dir, rel_root)
        mkpath(target_root)
        for file in files
            cp(joinpath(root, file), joinpath(target_root, file); force = true)
        end
    end
    return String(target_dir)
end

function browser_open_command(url::AbstractString)
    if Sys.isapple()
        return `open $url`
    elseif Sys.iswindows()
        return `cmd /c start "" $url`
    else
        return `xdg-open $url`
    end
end

function open_dashboard_browser(url::AbstractString)
    try
        run(ignorestatus(browser_open_command(url)))
    catch err
        @warn "Could not open dashboard in the default browser" url exception = (err, catch_backtrace())
    end
end

function dashboard_source_path(source::Union{Nothing, AbstractString})
    source === nothing && return nothing
    return basename(String(source))
end

function ensure_sized_dashboard_aircraft(ac::aircraft, label::AbstractString)
    ac.is_sized[1] || throw(ArgumentError("$label must be sized before exporting dashboard data"))
    return ac
end

function chord_eta(ac::aircraft)
    wing = ac.wing
    bo = wing.layout.root_span
    bs = wing.layout.break_span
    b = wing.layout.span
    co = wing.layout.root_chord
    cs = co * wing.inboard.λ
    ct = co * wing.outboard.λ
    eta = [0.0, safe_ratio(bo, b), safe_ratio(bs, b), 1.0]
    chord = [safe_float(co), safe_float(co), safe_float(cs), safe_float(ct)]
    return eta, chord
end

function spanwise_cl(ac::aircraft, para)
    wing = ac.wing
    eta = [0.0, safe_ratio(wing.layout.root_span, wing.layout.span), safe_ratio(wing.layout.break_span, wing.layout.span), 1.0]
    cl_by_point = Vector{Vector{Float64}}(undef, iptotal)
    for ip in 1:iptotal
        clo = safe_float(para[iaclpo, ip])
        cls = safe_float(para[iaclps, ip])
        clt = safe_float(para[iaclpt, ip])
        cl_by_point[ip] = [clo, clo, cls, clt]
    end
    return eta, cl_by_point
end

function extract_structural_weights(ac::aircraft)
    parg = ac.parg
    fuse = ac.fuselage
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    lg = ac.landing_gear

    wbox = wing.inboard.webs.weight.W + wing.inboard.caps.weight.W
    wwing = wbox * (1 + wing.weight_frac_flap + wing.weight_frac_slat +
                    wing.weight_frac_ailerons + wing.weight_frac_leading_trailing_edge +
                    wing.weight_frac_ribs + wing.weight_frac_spoilers +
                    wing.weight_frac_attachments)
    wadd = parg[igWMTO] * fuse.HPE_sys.W + lg.nose_gear.weight.W + lg.main_gear.weight.W

    labels = ["Fuselage", "Wing", "H-Tail", "V-Tail", "Installed engines", "Fuel Tank", "Added"]
    values_t = tonnes_from_force.([
        fuse.weight, wwing, htail.weight, vtail.weight,
        parg[igWeng], parg[igWftank], wadd,
    ])
    return labels, values_t
end

function wing_load_payload(ac::aircraft)
    wing = ac.wing
    etao = safe_float(wing.layout.ηo)
    etas = safe_float(wing.layout.ηs)
    b = safe_float(wing.layout.span)

    if wing.has_strut || b <= 0.0
        eta_range = vcat(collect(range(etao, etas, length = 20)), collect(range(etas, 1.0, length = 20)))
        shear = fill(safe_float(wing.inboard.max_shear_load), length(eta_range))
        moment = fill(safe_float(wing.inboard.moment), length(eta_range))
        return (
            eta = [safe_float(value) for value in eta_range],
            shear_N = safe_array(shear),
            moment_Nm = safe_array(moment),
            markers_eta = (
                wing_box_end = safe_float(etao),
                planform_break = safe_float(etas),
            ),
        )
    end

    para = view(ac.para, :, :, 1)
    parg = ac.parg
    ip = ipcruise1

    γt = safe_float(wing.outboard.λ * para[iarclt, ip])
    γs = safe_float(wing.inboard.λ * para[iarcls, ip])
    Nlift = safe_float(parg[igNlift])
    Lhtail = safe_float(parg[igWMTO] * ac.htail.CL_CLmax * ac.htail.layout.S / max(wing.layout.S, eps(Float64)))
    po = safe_float(wing_loading(wing, para[iarclt, ip], para[iarcls, ip], Nlift, parg[igWMTO], Lhtail))

    dLt = safe_float(wing.tip_lift_loss * po * wing.layout.root_chord * γt * wing.outboard.λ)
    half_span = 0.5 * b

    function q_in(eta)
        return 0.5 * po * b * (1.0 + (γs - 1.0) * (eta - etao) / max(etas - etao, eps(Float64)))
    end

    function q_out(eta)
        return 0.5 * po * b * (γs + (γt - γs) * (eta - etas) / max(1.0 - etas, eps(Float64)))
    end

    engine_step_N = 0.0
    if ac.options.opt_engine_location == EngineLocation.Wing && Int(round(parg[igneng])) > 0
        eng_per_side = safe_float(parg[igWeng] / max(parg[igneng], 1.0))
        engine_step_N = safe_float(Nlift * eng_per_side)
    end

    λs = safe_float(wing.inboard.λ)
    λt = safe_float(wing.outboard.λ)

    function chord_in(eta)
        return 1.0 + (λs - 1.0) * (eta - etao) / max(etas - etao, eps(Float64))
    end

    function chord_out(eta)
        return λs + (λt - λs) * (eta - etas) / max(1.0 - etas, eps(Float64))
    end

    i_in = max((etas - etao) * (1.0 + λs + λs^2) / 3.0, eps(Float64))
    i_out = max((1.0 - etas) * (λs^2 + λs * λt + λt^2) / 3.0, eps(Float64))

    w_in_total = safe_float(Nlift * wing.inboard.weight)
    w_out_total = safe_float(Nlift * wing.outboard.weight)

    function q_weight_in(eta)
        return w_in_total * chord_in(eta)^2 / i_in
    end

    function q_weight_out(eta)
        return w_out_total * chord_out(eta)^2 / i_out
    end

    function segment_integral(f, a, b; n = 60)
        b <= a && return 0.0
        xs = collect(range(a, b, length = n))
        ys = f.(xs)
        total = 0.0
        for i in 1:(length(xs) - 1)
            total += 0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i])
        end
        return safe_float(total)
    end

    function segment_moment(f, station, a, b; n = 60)
        b <= a && return 0.0
        xs = collect(range(a, b, length = n))
        ys = [f(x) * half_span * (x - station) for x in xs]
        total = 0.0
        for i in 1:(length(xs) - 1)
            total += 0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i])
        end
        return safe_float(total)
    end

    function shear_at(eta)
        force = 0.0
        force += segment_integral(q_out, max(eta, etas), 1.0)
        force -= segment_integral(q_weight_out, max(eta, etas), 1.0)
        if eta < etas
            force += segment_integral(q_in, max(eta, etao), etas)
            force -= segment_integral(q_weight_in, max(eta, etao), etas)
        end
        eta < 1.0 && (force += dLt)
        eta <= etas && (force -= engine_step_N)
        return safe_float(force)
    end

    function moment_at(eta)
        bending = 0.0
        bending += segment_moment(q_out, eta, max(eta, etas), 1.0)
        bending -= segment_moment(q_weight_out, eta, max(eta, etas), 1.0)
        if eta < etas
            bending += segment_moment(q_in, eta, max(eta, etao), etas)
            bending -= segment_moment(q_weight_in, eta, max(eta, etao), etas)
        end
        eta < 1.0 && (bending += dLt * half_span * (1.0 - eta))
        eta <= etas && (bending -= engine_step_N * half_span * (etas - eta))
        return safe_float(bending)
    end

    eta_inboard = collect(range(etao, etas, length = 28))
    eta_outboard = collect(range(etas, 1.0, length = 28))
    eta_range = copy(eta_inboard)
    shear = [shear_at(eta) for eta in eta_inboard]
    moment = [moment_at(eta) for eta in eta_inboard]
    if !isempty(shear)
        shear[end] = safe_float(wing.outboard.max_shear_load - engine_step_N)
        moment[end] = safe_float(wing.outboard.moment)

        Δη = max(etas - etao, eps(Float64))
        root_shear_target = safe_float(wing.inboard.max_shear_load)
        root_moment_target = safe_float(wing.inboard.moment)
        ds0 = root_shear_target - shear[1]
        dm0 = root_moment_target - moment[1]
        k = 2.0 * dm0 / max(b * Δη, eps(Float64))
        v = 3.0 * ds0 - 6.0 * k
        u = ds0 - v

        for i in eachindex(eta_inboard)
            t = clamp((etas - eta_inboard[i]) / Δη, 0.0, 1.0)
            shear[i] += u * t + v * t^2
            moment[i] += 0.5 * b * Δη * (u * t^2 / 2.0 + v * t^3 / 3.0)
        end
    end

    if !isempty(eta_outboard)
        push!(eta_range, etas)
        push!(shear, safe_float(wing.outboard.max_shear_load))
        push!(moment, safe_float(wing.outboard.moment))
        for eta in eta_outboard[2:end]
            push!(eta_range, eta)
            push!(shear, shear_at(eta))
            push!(moment, moment_at(eta))
        end
    end

    return (
        eta = [safe_float(value) for value in eta_range],
        shear_N = safe_array(shear),
        moment_Nm = safe_array(moment),
        markers_eta = (
            wing_box_end = safe_float(etao),
            planform_break = safe_float(etas),
        ),
    )
end

function engine_weight_breakdown_payload(ac::aircraft)
    parg = ac.parg
    total_t = tonnes_from_force(parg[igWeng])
    hx_t = tonnes_from_force(parg[igWHXs])
    bare_no_hx_t = tonnes_from_force(parg[igWebare] - parg[igWHXs])
    nacelle_t = tonnes_from_force(parg[igWnace])
    support_t = tonnes_from_force(parg[igWeng] - parg[igWebare] - parg[igWnace])

    return (
        labels = ["Bare engine", "Heat exchangers", "Nacelle", "Add'l. + pylon"],
        values_t = [bare_no_hx_t, hx_t, nacelle_t, support_t],
        total_t = total_t,
    )
end

function compressor_speed_lines(map::engine.CompressorMap, piD::Real)
    abs(piD - 1.0) < 1.0e-9 && return Any[]

    lines = Any[]
    for i in eachindex(map.NcMap)
        wc_row = map.WcMap[i, :]
        pr_row = map.PRMap[i, :]
        mb_frac = [safe_float(value / map.defaults.Wc) for value in wc_row]
        pr_abs = @. 1.0 + (pr_row - 1.0) * (piD - 1.0) / (map.defaults.PR - 1.0)
        pr_frac = [safe_float(value / piD) for value in pr_abs]
        push!(lines, (
            x = mb_frac,
            y = pr_frac,
            speed_frac = safe_float(map.NcMap[i] / map.defaults.Nc),
        ))
    end
    return lines
end

function compressor_efficiency_grid(map::engine.CompressorMap, piD::Real, mbD::Real, NbD::Real, ep0::Real; n = 34)
    mb_min = min(0.25, minimum(map.WcMap) / map.defaults.Wc)
    mb_max = max(1.50, maximum(map.WcMap) / map.defaults.Wc)
    mb_grid = collect(range(mb_min, mb_max, length = n))
    pr_abs_map = @. 1.0 + (map.PRMap - 1.0) * (piD - 1.0) / (map.defaults.PR - 1.0)
    pr_min = minimum(pr_abs_map) / piD
    pr_max = max(1.50, maximum(pr_abs_map) / piD)
    pr_grid = collect(range(pr_min, pr_max, length = n))
    eff = Matrix{Float64}(undef, n, n)

    for (j, pr_frac) in enumerate(pr_grid)
        for (i, mb_frac) in enumerate(mb_grid)
            try
                _, epol, _, _, _, _, _, _ = engine.calculate_compressor_speed_and_efficiency(
                    map, pr_frac * piD, mb_frac * mbD, piD, mbD, NbD, ep0)
                eff[j, i] = epol
            catch
                eff[j, i] = NaN
            end
        end
    end

    return (
        x = [safe_float(value) for value in mb_grid],
        y = [safe_float(value) for value in pr_grid],
        z = safe_nested(eff),
    )
end

function compressor_map_payload(map::engine.CompressorMap, pare;
        label::AbstractString,
        pr_idx::Int,
        mb_idx::Int,
        nb_idx::Int,
        pr_design_idx::Int,
        mb_design_idx::Int,
        nb_design_idx::Int,
        epol0_idx::Int,
        eff_idx::Int)
    piD = safe_float(pare[pr_design_idx, ipcruise1])
    mbD = safe_float(pare[mb_design_idx, ipcruise1])
    NbD = safe_float(pare[nb_design_idx, ipcruise1])
    ep0 = safe_float(pare[epol0_idx, ipcruise1])

    (piD <= 1.0 || mbD <= 0.0 || NbD <= 0.0 || ep0 <= 0.0) && return nothing

    return (
        label = String(label),
        design = (
            pressure_ratio = piD,
            corrected_mass_flow = mbD,
            corrected_speed = NbD,
            polytropic_efficiency = ep0,
        ),
        speed_lines = compressor_speed_lines(map, piD),
        efficiency_grid = compressor_efficiency_grid(map, piD, mbD, NbD, ep0),
        operating_line = (
            x = [safe_ratio(pare[mb_idx, ip], mbD) for ip in 1:iptotal],
            y = [safe_ratio(pare[pr_idx, ip], piD) for ip in 1:iptotal],
            speed_frac = [safe_ratio(pare[nb_idx, ip], NbD) for ip in 1:iptotal],
            eff = [safe_float(pare[eff_idx, ip]) for ip in 1:iptotal],
        ),
    )
end

function wing_polygon_points(ac::aircraft)
    wing = ac.wing
    co = wing.layout.root_chord
    cs = co * wing.inboard.λ
    ct = co * wing.outboard.λ
    bo = wing.layout.root_span
    bs = wing.layout.break_span
    b = wing.layout.span
    dx = wing.layout.box_x
    tanl = tan(wing.layout.sweep * pi / 180.0)
    xax = 0.40

    xs = tanl * (bs - bo) / 2.0
    xt = tanl * (b - bo) / 2.0

    xw = [
        co * (-xax) + dx,
        xs + cs * (-xax) + dx,
        xt + ct * (-xax) + dx,
        xt + ct * (1 - xax) + dx,
        xs + cs * (1 - xax) + dx,
        co * (1 - xax) + dx,
    ]
    yw = [bo / 2, bs / 2, b / 2, b / 2, bs / 2, bo / 2]

    xs_all = vcat(xw, reverse(xw))
    ys_all = vcat(yw, reverse(-yw))
    return [[safe_float(x), safe_float(y)] for (x, y) in zip(xs_all, ys_all)]
end

function fuselage_polygon_points(ac::aircraft)
    fuselage = ac.fuselage
    radius = fuselage.layout.radius
    bubble_offset = fuselage.layout.bubble_center_y_offset
    nose_radius = max(fuselage.layout.nose_radius, 0.5)
    tail_radius = max(fuselage.layout.tail_radius, 0.5)
    x0 = fuselage.layout.x_nose
    x1 = fuselage.layout.x_start_cylinder
    x2 = fuselage.layout.x_end_cylinder
    xe = fuselage.layout.x_end
    half_width = radius + bubble_offset
    dy = startswith(string(fuselage.layout.opt_tapers_to), "p") ? -half_width : -0.2 * half_width
    npts = 24

    xf = zeros(npts + npts + 1)
    yf = zeros(npts + npts + 1)

    for i in 1:npts
        t = (i - 1) / (npts - 1)
        fx = cos(0.5 * pi * t)
        xf[i] = x1 + (x0 - x1) * fx
        yf[i] = half_width * (1 - fx^nose_radius)^(1 / nose_radius)
    end

    for i in 1:npts
        t = (i - 1) / (npts - 1)
        xf[npts + i] = x2 + (xe - x2) * t
        yf[npts + i] = half_width + dy * t^tail_radius
    end

    xf[end] = xf[end - 1]
    yf[end] = 0.0

    xs = vcat(xf, reverse(xf))
    ys = vcat(yf, reverse(-yf))
    return [[safe_float(x), safe_float(y)] for (x, y) in zip(xs, ys)]
end

function htail_polygon_points(ac::aircraft)
    htail = ac.htail
    fuselage = ac.fuselage
    sh = htail.layout.S
    arh = htail.layout.AR
    lambdah = htail.outboard.λ
    boh = htail.layout.root_span
    sweeph = htail.layout.sweep
    dx = htail.layout.box_x
    bh = sqrt(sh * arh)
    coh = sh / (boh + (bh - boh) * 0.5 * (1 + lambdah))
    cth = coh * lambdah
    tanl = tan(sweeph * pi / 180.0)
    xax = 0.40

    xole = coh * (-xax) + dx
    xote = coh * (1 - xax) + dx
    xtle = cth * (-xax) + dx + 0.5 * (bh - boh) * tanl
    xtte = cth * (1 - xax) + dx + 0.5 * (bh - boh) * tanl

    if startswith(string(fuselage.layout.opt_tapers_to), "p")
        xcle = xole
        xcte = xote
        ycle = 0.5 * boh
        ycte = 0.5 * boh
    else
        xcle = coh * (-xax) + dx - 0.5 * boh * tanl
        xcte = coh * (1 - xax) + dx - 0.5 * boh * tanl
        ycle = 0.0
        ycte = 0.0
    end

    xh = [xcle, xole, xtle, xtte, xote, xcte]
    yh = [ycle, 0.5 * boh, 0.5 * bh, 0.5 * bh, 0.5 * boh, ycte]
    xs = vcat(xh, reverse(xh))
    ys = vcat(yh, reverse(-yh))
    return [[safe_float(x), safe_float(y)] for (x, y) in zip(xs, ys)]
end

function engine_rectangles(ac::aircraft)
    parg = ac.parg
    wing = ac.wing
    diameter = parg[igdfan]
    neng = Int(parg[igneng])
    nacelle_length = parg[iglnace]
    span = wing.layout.span
    bo = wing.layout.root_span
    co = wing.layout.root_chord
    lambdas = wing.inboard.λ
    lambdat = wing.outboard.λ
    etas = wing.layout.ηs
    tanl = tan(wing.layout.sweep * pi / 180.0)
    dx = wing.layout.box_x

    dy = 2 * diameter
    yi = if neng == 2
        [etas * span / 2]
    else
        collect(range(bo / 2 + dy, stop = span / 2 * 3 / 4, length = max(1, Int(neng / 2))))
    end

    rects = Vector{Vector{Vector{Float64}}}()
    for y in yi
        eta = y / (span / 2)
        etao = bo / span
        ci = if eta <= etas
            co * (1 + (lambdas - 1) * (eta - etao) / (etas - etao))
        else
            co * (lambdas + (lambdat - lambdas) * (eta - etas) / (1 - etas))
        end
        xi = tanl * (y - bo / 2) - 0.4 * ci + dx - 1.0
        for sign in (-1.0, 1.0)
            yc = sign * y
            push!(rects, [
                [safe_float(xi), safe_float(yc - diameter / 2)],
                [safe_float(xi + nacelle_length), safe_float(yc - diameter / 2)],
                [safe_float(xi + nacelle_length), safe_float(yc + diameter / 2)],
                [safe_float(xi), safe_float(yc + diameter / 2)],
            ])
        end
    end
    return rects
end

function component_bounds(points)
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    return minimum(xs), maximum(xs), minimum(ys), maximum(ys)
end

function bubble_center_positions(layout)
    layout.n_webs == 0 && return [0.0]
    return [safe_float(multiplier * layout.bubble_center_y_offset) for multiplier in (-layout.n_webs):2:layout.n_webs]
end

function fuselage_cross_section_points(fuselage; npts = 220)
    layout = fuselage.layout
    radius = safe_float(layout.radius)
    downward_shift = safe_float(layout.bubble_lower_downward_shift)
    centers = bubble_center_positions(layout)

    xmin = minimum(centers) - radius
    xmax = maximum(centers) + radius
    xs = collect(range(xmin, xmax, length = npts))

    upper = [
        maximum(sqrt(max(radius^2 - (x - center)^2, 0.0)) for center in centers)
        for x in xs
    ]
    lower_xs = reverse(xs)
    lower = [
        -downward_shift - maximum(sqrt(max(radius^2 - (x - center)^2, 0.0)) for center in centers)
        for x in lower_xs
    ]

    top_points = [[safe_float(x), safe_float(y)] for (x, y) in zip(xs, upper)]
    bottom_points = [[safe_float(x), safe_float(y)] for (x, y) in zip(lower_xs, lower)]
    return vcat(top_points, bottom_points)
end

function fuselage_cabin_width(ac::aircraft)
    ac.options.is_doubledecker && return nothing

    fuselage = ac.fuselage
    try
        theta = find_floor_angles(false, fuselage.layout.radius, fuselage.layout.bubble_lower_downward_shift;
            h_seat = fuselage.cabin.seat_height)
        width = find_cabin_width(
            fuselage.layout.radius,
            fuselage.layout.bubble_center_y_offset,
            fuselage.layout.n_webs,
            theta,
            fuselage.cabin.seat_height,
        )
        return safe_float(width)
    catch
        return nothing
    end
end

function trefftz_payload(ac::aircraft)
    try
        trefftz_config = ac.options.trefftz_config
        i_w1 = aerodynamics.i_first_wing(trefftz_config)
        i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1
        i_t1 = aerodynamics.i_first_tail(trefftz_config)
        i_t2 = aerodynamics.i_last_tail(trefftz_config)

        circulation_by_point = Vector{Vector{Float64}}(undef, iptotal)
        wake_circulation_by_point = Vector{Vector{Float64}}(undef, iptotal)
        induced_velocity_mps = Vector{Vector{Float64}}(undef, iptotal)
        drag_density_proxy = Vector{Vector{Float64}}(undef, iptotal)
        wing_span_y_m = Float64[]
        wake_span_y_m = Float64[]

        for ip in 1:iptotal
            para_view = @view ac.para[:, ip, 1]
            aerodynamics.induced_drag!(para_view, ac, trefftz_config)

            ws = ac.wake_system
            ws === nothing && return nothing

            ycp = aerodynamics.ctrl_ys(ws)
            ys = aerodynamics.field_ys(ws)
            gc = aerodynamics.TREFFTZ_GEOM.gc[1:length(ycp)]

            if isempty(wing_span_y_m)
                wing_span_y_m = [safe_float(y) for y in ycp[i_w1:i_w2]]
                wake_span_y_m = [safe_float(y) for y in ys[i_w1:(i_w2 + 1)]]
            end

            gw = zeros(length(ys))
            ifrst = [i_w1, i_t1]
            ilast = [i_w2 + 1, i_t2]
            aerodynamics.calculate_wake_circulation!(gw, gc, ifrst, ilast, 2)

            vnc = ws.influence_matrix * gw
            vnc .*= ac.wing.layout.span / (2.0π)

            circulation_by_point[ip] = [safe_float(value) for value in gc[i_w1:i_w2]]
            wake_circulation_by_point[ip] = [safe_float(value) for value in gw[i_w1:(i_w2 + 1)]]
            induced_velocity_mps[ip] = [safe_float(value) for value in vnc[i_w1:i_w2]]
            drag_density_proxy[ip] = [safe_float(-gc[i] * vnc[i]) for i in i_w1:i_w2]
        end

        isempty(wing_span_y_m) && return nothing

        return (;
            wing_span_y_m = wing_span_y_m,
            wake_span_y_m = wake_span_y_m,
            circulation_by_point = circulation_by_point,
            wake_circulation_by_point = wake_circulation_by_point,
            induced_velocity_mps = induced_velocity_mps,
            drag_density_proxy = drag_density_proxy,
            span_markers_m = (;
                root = safe_float(ac.wing.layout.root_span / 2),
                span_break = safe_float(ac.wing.layout.break_span / 2),
                tip = safe_float(ac.wing.layout.span / 2),
            ),
        )
    catch err
        @debug "Skipping Trefftz dashboard payload" aircraft=ac.name exception=(err, catch_backtrace())
        return nothing
    end
end

function geometry_payload(ac::aircraft, para)
    wing = ac.wing
    fuselage = ac.fuselage
    htail = ac.htail
    parg = ac.parg
    cabin_width_m = fuselage_cabin_width(ac)

    wing_points = wing_polygon_points(ac)
    fuselage_points = fuselage_polygon_points(ac)
    htail_points = htail_polygon_points(ac)
    engine_points = engine_rectangles(ac)
    cross_section_points = fuselage_cross_section_points(fuselage)
    bubble_centers = bubble_center_positions(fuselage.layout)

    all_points = vcat(wing_points, fuselage_points, htail_points, reduce(vcat, engine_points; init = Vector{Vector{Float64}}()))
    xs = [p[1] for p in all_points]
    ys = [p[2] for p in all_points]
    xpad = max(2.0, 0.05 * (maximum(xs) - minimum(xs)))
    ypad = max(2.0, 0.10 * (maximum(ys) - minimum(ys)))
    view_box = (
        x = safe_float(minimum(xs) - xpad),
        y = safe_float(minimum(ys) - ypad),
        width = safe_float(maximum(xs) - minimum(xs) + 2 * xpad),
        height = safe_float(maximum(ys) - minimum(ys) + 2 * ypad),
    )

    wing_lod = safe_ratio(para[iaCL, ipcruise1], para[iaCD, ipcruise1])
    wing_tip = component_bounds(wing_points)
    fuselage_bounds = component_bounds(fuselage_points)
    tail_bounds = component_bounds(htail_points)
    cross_section_bounds = component_bounds(cross_section_points)
    fuselage_width = safe_float(cross_section_bounds[2] - cross_section_bounds[1])
    fuselage_height = safe_float(cross_section_bounds[4] - cross_section_bounds[3])
    section_label = fuselage.layout.n_webs > 0 ? "Multi bubble" : "Single bubble"
    cabin_width_text = cabin_width_m === nothing ? "n/a" : @sprintf("%.2f m", cabin_width_m)

    components = Any[
        (
            kind = "wing",
            points = wing_points,
            label = "Wing",
            tooltip = @sprintf("Wing\nSpan: %.1f m\nArea: %.1f m^2\nAR: %.2f\nSweep: %.1f deg\nCruise L/D: %.2f",
                wing.layout.span, wing.layout.S, wing.layout.AR, wing.layout.sweep, wing_lod),
        ),
        (
            kind = "fuselage",
            points = fuselage_points,
            label = "Fuselage",
            tooltip = @sprintf("Fuselage\nSection: %s\nLength: %.1f m\nWidth: %.2f m\nCabin width: %s\nWeight: %.2f t",
                section_label, fuselage.layout.x_end - fuselage.layout.x_nose, fuselage_width, cabin_width_text, tonnes_from_force(fuselage.weight)),
        ),
        (
            kind = "htail",
            points = htail_points,
            label = "H-Tail",
            tooltip = @sprintf("Horizontal tail\nArea: %.1f m^2\nAR: %.2f\nSweep: %.1f deg",
                htail.layout.S, htail.layout.AR, htail.layout.sweep),
        ),
    ]

    for rect in engine_points
        push!(components, (
            kind = "engine",
            points = rect,
            label = "Engine",
            tooltip = @sprintf("Engine nacelle\nFan diameter: %.2f m\nNacelle length: %.2f m\nEngines: %d",
                parg[igdfan], parg[iglnace], Int(parg[igneng])),
        ))
    end

    chord_eta_values, chord_values = chord_eta(ac)

    return (
        wing = (
            span_m = safe_float(wing.layout.span),
            area_m2 = safe_float(wing.layout.S),
            root_chord_m = safe_float(wing.layout.root_chord),
            aspect_ratio = safe_float(wing.layout.AR),
            sweep_deg = safe_float(wing.layout.sweep),
            chord = (
                eta = chord_eta_values,
                chord_m = chord_values,
            ),
        ),
        fuselage = (
            radius_m = safe_float(fuselage.layout.radius),
            length_m = safe_float(fuselage.layout.x_end - fuselage.layout.x_nose),
            width_m = fuselage_width,
            height_m = fuselage_height,
            cabin_width_m = cabin_width_m,
            cross_section = (
                type = fuselage.layout.n_webs > 0 ? "multi_bubble" : "single_bubble",
                points = cross_section_points,
                bubble_centers_y_m = bubble_centers,
            ),
        ),
        svg = (
            view_box = view_box,
            components = components,
            focus_points = (
                wing = (x = safe_float((wing_tip[1] + wing_tip[2]) / 2), y = safe_float(wing_tip[4] * 0.72)),
                fuselage = (x = safe_float((fuselage_bounds[1] + fuselage_bounds[2]) / 2), y = 0.0),
                htail = (x = safe_float((tail_bounds[1] + tail_bounds[2]) / 2), y = safe_float(tail_bounds[4] * 0.82)),
            ),
        ),
    )
end

function summary_payload(ac::aircraft, para, pare)
    parg = ac.parg
    parm = ac.parm
    mtow_t = tonnes_from_force(parg[igWMTO])
    fuel_t = tonnes_from_force(parg[igWfuel])
    payload_t = tonnes_from_force(parg[igWpay])
    oew_t = mtow_t - fuel_t - payload_t
    cruise_lod = safe_ratio(para[iaCL, ipcruise1], para[iaCD, ipcruise1])

    return (
        pfei_kj_kg_km = safe_float(parm[imPFEI, 1]),
        range_km = safe_float(parm[imRange, 1] / 1e3),
        cruise_lod = cruise_lod,
        cruise_mach = safe_float(para[iaMach, ipcruise1]),
        cruise_tsfc_mg_ns = kg_per_ns_to_mg_per_ns(pare[ieTSFC, ipcruise1]),
        mtow_t = mtow_t,
        oew_t = safe_float(oew_t),
        fuel_t = fuel_t,
        payload_t = payload_t,
        wing_span_m = safe_float(ac.wing.layout.span),
        wing_area_m2 = safe_float(ac.wing.layout.S),
        wing_ar = safe_float(ac.wing.layout.AR),
        wing_sweep_deg = safe_float(ac.wing.layout.sweep),
        fan_diameter_m = safe_float(parg[igdfan]),
        cruise_bpr = safe_float(pare[ieBPR, ipcruise1]),
        cruise_opr = safe_float(pare[ieOPR, ipcruise1]),
        cruise_tt4_k = safe_float(pare[ieTt4, ipcruise1]),
        kpis = Any[
            (key = "pfei_kj_kg_km", label = raw"\(\mathrm{PFEI}\) [kJ/kg/km]", value = safe_float(parm[imPFEI, 1]), digits = 3, lower_better = true),
            (key = "range_km", label = raw"\(R_{des}\) [km]", value = safe_float(parm[imRange, 1] / 1e3), digits = 0, lower_better = false),
            (key = "cruise_lod", label = raw"Cruise \(L/D\)", value = cruise_lod, digits = 2, lower_better = false),
            (key = "cruise_mach", label = raw"Cruise \(M\)", value = safe_float(para[iaMach, ipcruise1]), digits = 3, lower_better = false),
            (key = "cruise_tsfc_mg_ns", label = raw"Cruise \(\mathrm{TSFC}\) [mg/N/s]", value = kg_per_ns_to_mg_per_ns(pare[ieTSFC, ipcruise1]), digits = 4, lower_better = true),
            (key = "mtow_t", label = raw"\(W_{MTO}\) [t]", value = mtow_t, digits = 1, lower_better = true),
            (key = "oew_t", label = raw"\(W_{OE}\) [t]", value = safe_float(oew_t), digits = 1, lower_better = true),
            (key = "fuel_t", label = raw"\(W_f\) [t]", value = fuel_t, digits = 1, lower_better = true),
            (key = "payload_t", label = raw"\(W_{pay}\) [t]", value = payload_t, digits = 1, lower_better = false),
            (key = "wing_span_m", label = raw"\(b\) [m]", value = safe_float(ac.wing.layout.span), digits = 2, lower_better = false),
            (key = "wing_ar", label = raw"Wing \(AR\)", value = safe_float(ac.wing.layout.AR), digits = 2, lower_better = false),
            (key = "wing_sweep_deg", label = raw"Wing \(\Lambda\) [deg]", value = safe_float(ac.wing.layout.sweep), digits = 1, lower_better = false),
            (key = "fan_diameter_m", label = raw"Fan \(D_f\) [m]", value = safe_float(parg[igdfan]), digits = 3, lower_better = false),
            (key = "cruise_bpr", label = raw"Cruise \(\mathrm{BPR}\)", value = safe_float(pare[ieBPR, ipcruise1]), digits = 2, lower_better = false),
            (key = "cruise_opr", label = raw"Cruise \(\mathrm{OPR}\)", value = safe_float(pare[ieOPR, ipcruise1]), digits = 1, lower_better = false),
            (key = "cruise_tt4_k", label = raw"Cruise \(T_{t4}\) [K]", value = safe_float(pare[ieTt4, ipcruise1]), digits = 0, lower_better = true),
        ],
    )
end

function mission_payload(para)
    return (
        labels = DASHBOARD_MISSION_LABELS,
        short_labels = DASHBOARD_MISSION_SHORT,
        range_km = safe_array(para[iaRange, :] ./ 1e3),
        time_hr = safe_array(para[iatime, :] ./ 3600),
        gamma_rad = safe_array(para[iagamV, :]),
        alt_km = safe_array(para[iaalt, :] ./ 1e3),
        mach = safe_array(para[iaMach, :]),
        lod = [safe_ratio(para[iaCL, ip], para[iaCD, ip]) for ip in 1:iptotal],
        roc_mps = safe_array(para[iaROC, :]),
        weight_fraction = safe_array(para[iafracW, :]),
        cd_e4 = safe_array(para[iaCD, :] .* 1e4),
    )
end

function aeroperf_cl_range(para)
    cls = Float64[safe_float(para[iaCL, ip]) for ip in 1:iptotal if para[iaCL, ip] isa Real && isfinite(para[iaCL, ip]) && para[iaCL, ip] > 0.05]
    if isempty(cls)
        return collect(range(0.15, 0.90; length = 16))
    end

    cl_min = clamp(minimum(cls) * 0.70, 0.12, 0.55)
    cl_max = clamp(maximum(cls) * 1.10, cl_min + 0.20, 0.90)
    return collect(range(cl_min, cl_max; length = 16))
end

function aeroperf_mach_values(para)
    cruise_mach = safe_float(para[iaMach, ipcruise1])
    candidates = [round(clamp(cruise_mach + ΔM, 0.30, 0.90); digits = 3) for ΔM in (-0.05, 0.0, 0.05)]
    return unique(candidates)
end

function aeroperf_sweep_payload(ac::aircraft, para)
    cl_range = aeroperf_cl_range(para)
    mach_values = aeroperf_mach_values(para)
    curves = Any[]

    for mach in mach_values
        try
            sweep = aeroperf_sweep(ac, cl_range; Mach = mach, ip = ipcruise1)
            field_values = (; (field.key => safe_array(getfield(sweep, field.key)) for field in DASHBOARD_AEROPERF_FIELDS)...)
            push!(curves, merge((
                mach = safe_float(sweep.Mach),
                CLs = safe_array(sweep.CLs),
            ), field_values))
        catch
            continue
        end
    end

    isempty(curves) && return nothing

    return (
        default_field = "LDs",
        fields = [(key = String(field.key), label = field.label) for field in DASHBOARD_AEROPERF_FIELDS],
        curves = curves,
    )
end

function aero_payload(ac::aircraft, para)
    area = ac.wing.layout.S
    chord_eta_values, chord_values = chord_eta(ac)
    cl_eta, cl_by_point = spanwise_cl(ac, para)
    drag_fracs = [
        [safe_ratio(para[idx, ip], para[iaCD, ip]) for ip in 1:iptotal]
        for idx in DASHBOARD_DRAG_INDICES
    ]

    return (
        drag_keys = DASHBOARD_DRAG_KEYS,
        drag_area_cruise_cm2 = [
            safe_float(para[idx, ipcruise1] * area * 1e4) for idx in DASHBOARD_DRAG_INDICES
        ],
        drag_fractions = drag_fracs,
        lod = [safe_ratio(para[iaCL, ip], para[iaCD, ip]) for ip in 1:iptotal],
        tail_volume = (
            htail = safe_float(ac.htail.volume),
            vtail = safe_float(ac.vtail.volume),
        ),
        wing = (
            eta = chord_eta_values,
            chord_m = chord_values,
            cl_eta = cl_eta,
            cl_by_point = cl_by_point,
        ),
        aeroperf_sweep = aeroperf_sweep_payload(ac, para),
        trefftz = trefftz_payload(ac),
    )
end

function engine_payload(ac::aircraft, pare)
    parg = ac.parg
    tt_by_point = Vector{Vector{Float64}}(undef, iptotal)
    pt_by_point = Vector{Vector{Float64}}(undef, iptotal)
    for ip in 1:iptotal
        tt_by_point[ip] = [safe_float(pare[idx, ip]) for idx in DASHBOARD_T_STATION_INDICES]
        pt_by_point[ip] = [safe_float(pare[idx, ip] / 1e3) for idx in DASHBOARD_P_STATION_INDICES]
    end

    return (
        temperature_stations = DASHBOARD_T_STATION_LABELS,
        pressure_stations = DASHBOARD_P_STATION_LABELS,
        temperature_k = tt_by_point,
        pressure_kpa = pt_by_point,
        fan_diameter_m = safe_float(parg[igdfan]),
        thrust_total_kN = safe_array(pare[ieFe, :] .* parg[igneng] ./ 1e3),
        thrust_per_engine_kN = safe_array(pare[ieFe, :] ./ 1e3),
        tsfc_mg_ns = safe_array(pare[ieTSFC, :] .* 1e6),
        bpr = safe_array(pare[ieBPR, :]),
        opr = safe_array(pare[ieOPR, :]),
        fpr = safe_array(pare[iepif, :]),
        fan_pr = safe_array(pare[iepif, :]),
        lpc_pr = safe_array(pare[iepilc, :]),
        hpc_pr = safe_array(pare[iepihc, :]),
        fan_eff = safe_array(pare[ieepf, :]),
        lpc_eff = safe_array(pare[ieeplc, :]),
        hpc_eff = safe_array(pare[ieephc, :]),
        fan_speed_frac = [safe_ratio(pare[ieNbf, ip], pare[ieNbfD, ipcruise1]) for ip in 1:iptotal],
        lpc_speed_frac = [safe_ratio(pare[ieNblc, ip], pare[ieNblcD, ipcruise1]) for ip in 1:iptotal],
        hpc_speed_frac = [safe_ratio(pare[ieNbhc, ip], pare[ieNbhcD, ipcruise1]) for ip in 1:iptotal],
        fan_mass_flow_frac = [safe_ratio(pare[iembf, ip], pare[iembfD, ipcruise1]) for ip in 1:iptotal],
        lpc_mass_flow_frac = [safe_ratio(pare[iemblc, ip], pare[iemblcD, ipcruise1]) for ip in 1:iptotal],
        hpc_mass_flow_frac = [safe_ratio(pare[iembhc, ip], pare[iembhcD, ipcruise1]) for ip in 1:iptotal],
        tt4_k = safe_array(pare[ieTt4, :]),
        weight_breakdown = engine_weight_breakdown_payload(ac),
        maps = (
            fan = compressor_map_payload(engine.FanMap, pare;
                label = "Fan",
                pr_idx = iepif,
                mb_idx = iembf,
                nb_idx = ieNbf,
                pr_design_idx = iepifD,
                mb_design_idx = iembfD,
                nb_design_idx = ieNbfD,
                epol0_idx = ieepolf,
                eff_idx = ieepf),
            lpc = compressor_map_payload(engine.LPCMap, pare;
                label = "LPC",
                pr_idx = iepilc,
                mb_idx = iemblc,
                nb_idx = ieNblc,
                pr_design_idx = iepilcD,
                mb_design_idx = iemblcD,
                nb_design_idx = ieNblcD,
                epol0_idx = ieepollc,
                eff_idx = ieeplc),
            hpc = compressor_map_payload(engine.HPCMap, pare;
                label = "HPC",
                pr_idx = iepihc,
                mb_idx = iembhc,
                nb_idx = ieNbhc,
                pr_design_idx = iepihcD,
                mb_design_idx = iembhcD,
                nb_design_idx = ieNbhcD,
                epol0_idx = ieepolhc,
                eff_idx = ieephc),
        ),
    )
end

function weights_payload(ac::aircraft)
    parg = ac.parg
    mtow_t = tonnes_from_force(parg[igWMTO])
    fuel_t = tonnes_from_force(parg[igWfuel])
    payload_t = tonnes_from_force(parg[igWpay])
    oew_t = mtow_t - fuel_t - payload_t
    labels, structure_values = extract_structural_weights(ac)

    return (
        mtow_t = mtow_t,
        oew_t = safe_float(oew_t),
        fuel_t = fuel_t,
        payload_t = payload_t,
        mtow_breakdown = (
            labels = ["OEW", "Fuel", "Payload"],
            values_t = [safe_float(oew_t), fuel_t, payload_t],
        ),
        structure_breakdown = (
            labels = labels,
            values_t = structure_values,
        ),
        wing_loads = wing_load_payload(ac),
    )
end

function dashboard_aircraft_payload(ac::aircraft, id::AbstractString; source::Union{Nothing, AbstractString} = nothing)
    ensure_sized_dashboard_aircraft(ac, id)

    para = view(ac.para, :, :, 1)
    pare = view(ac.pare, :, :, 1)

    return (
        id = String(id),
        name = String(ac.name),
        source = dashboard_source_path(source),
        summary = summary_payload(ac, para, pare),
        geometry = geometry_payload(ac, para),
        mission = mission_payload(para),
        weights = weights_payload(ac),
        aero = aero_payload(ac, para),
        engine = engine_payload(ac, pare),
        raw = raw_payload(ac),
    )
end

function dashboard_comparison(aircraft_payloads)
    length(aircraft_payloads) == 2 || return nothing
    left = aircraft_payloads[1]
    right = aircraft_payloads[2]
    rows = Any[]

    for (metric_left, metric_right) in zip(left.summary.kpis, right.summary.kpis)
        delta_abs = safe_float(metric_right.value - metric_left.value)
        delta_pct = abs(metric_left.value) < eps(Float64) ? nothing : safe_float(100 * delta_abs / abs(metric_left.value))
        push!(rows, (
            key = metric_left.key,
            label = metric_left.label,
            left = metric_left.value,
            right = metric_right.value,
            delta_abs = delta_abs,
            delta_pct = delta_pct,
            digits = metric_left.digits,
            lower_better = metric_left.lower_better,
        ))
    end

    return (
        left_name = left.name,
        right_name = right.name,
        summary = rows,
    )
end

function dashboard_payload(ac1::aircraft, ac2::Union{Nothing, aircraft} = nothing;
        source1::Union{Nothing, AbstractString} = nothing,
        source2::Union{Nothing, AbstractString} = nothing)
    aircraft_payloads = ac2 === nothing ?
        [dashboard_aircraft_payload(ac1, "ac1"; source = source1)] :
        [dashboard_aircraft_payload(ac1, "ac1"; source = source1), dashboard_aircraft_payload(ac2, "ac2"; source = source2)]

    return (
        schema_version = DASHBOARD_SCHEMA_VERSION,
        meta = (
            title = "TASOPT Diagnostic Dashboard",
            generated_at = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            default_point = ipcruise1,
            mission_labels = DASHBOARD_MISSION_LABELS,
            mission_short = DASHBOARD_MISSION_SHORT,
            plotly_src = DASHBOARD_PLOTLY_SRC,
        ),
        aircraft = aircraft_payloads,
        comparison = dashboard_comparison(aircraft_payloads),
    )
end

function load_dashboard_aircraft(ac_or_path::aircraft)
    return ac_or_path, nothing
end

function load_dashboard_aircraft(ac_or_path::AbstractString)
    ac = quickload_aircraft(String(ac_or_path))
    return ac, String(ac_or_path)
end

"""
    export_dashboard_data(ac1, ac2=nothing; filepath="dashboard.json")

Export chart-ready and SVG-ready TASOPT dashboard data to a single JSON file.
`ac1` and `ac2` may be either sized `aircraft` objects or paths to `.jld2`
files created with [`quicksave_aircraft`](@ref).
"""
function export_dashboard_data(ac1, ac2 = nothing; filepath::AbstractString = "dashboard.json")
    ac_left, source_left = load_dashboard_aircraft(ac1)
    ac_right, source_right = ac2 === nothing ? (nothing, nothing) : load_dashboard_aircraft(ac2)
    payload = dashboard_payload(ac_left, ac_right; source1 = source_left, source2 = source_right)
    json = _to_json(payload)
    write(filepath, json)
    return String(filepath)
end

function export_dashboard_data(ac1::aircraft, filepath::AbstractString)
    return export_dashboard_data(ac1; filepath = filepath)
end

function export_dashboard_data(ac1, ac2, filepath::AbstractString)
    return export_dashboard_data(ac1, ac2; filepath = filepath)
end

function injected_index_html(json::AbstractString)
    dist_dir, _ = ensure_dashboard_dist()
    dist_index = read(joinpath(dist_dir, "index.html"), String)
    data_tag = "<script id=\"taesopt-dashboard-data\" type=\"application/json\">" * json * "</script>\n"
    occursin("</body>", dist_index) || throw(ArgumentError(
        "Built dashboard HTML does not contain a closing </body> tag and cannot be injected."))
    return replace(dist_index, "</body>" => data_tag * "</body>")
end

"""
    export_dashboard_bundle(ac1, ac2=nothing; outdir="dashboard_bundle")

Copy the committed static dashboard assets into `outdir`, write `dashboard.json`,
and generate an `index.html` with the dashboard data embedded for direct browser
opening. This keeps the TypeScript/JS app static while letting Julia inject the
current aircraft snapshot.
"""
function export_dashboard_bundle(ac1, ac2 = nothing; outdir::AbstractString = "dashboard_bundle")
    ac_left, source_left = load_dashboard_aircraft(ac1)
    ac_right, source_right = ac2 === nothing ? (nothing, nothing) : load_dashboard_aircraft(ac2)
    payload = dashboard_payload(ac_left, ac_right; source1 = source_left, source2 = source_right)
    json = _to_json(payload)

    copy_dashboard_dist!(outdir)

    write(joinpath(outdir, "dashboard.json"), json)
    write(joinpath(outdir, "index.html"), injected_index_html(json))
    return joinpath(String(outdir), "index.html")
end

function export_dashboard_bundle(ac1::aircraft, outdir::AbstractString)
    return export_dashboard_bundle(ac1; outdir = outdir)
end

function export_dashboard_bundle(ac1, ac2, outdir::AbstractString)
    return export_dashboard_bundle(ac1, ac2; outdir = outdir)
end

"""
    start_dashboard(; open_browser=true)

Open the committed static dashboard shell from `dashboard/dist/index.html`.
The dashboard starts with no aircraft loaded; the user selects a dashboard JSON
file from the browser file picker.
"""
function start_dashboard(; open_browser::Bool = true)
    dist_dir, _ = ensure_dashboard_dist()
    index_path = joinpath(dist_dir, "index.html")
    open_browser && open_dashboard_browser(index_path)
    return index_path
end
