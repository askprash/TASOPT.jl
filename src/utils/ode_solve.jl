"""
    ode_solve(f, y0, tspan; tstops = Float64[], reltol = 1e-6, abstol = 1e-6, maxiters = 100_000)

Integrates the initial-value problem `dy/dt = f(y, t)`, `y(tspan[1]) = y0`, from `tspan[1]` to `tspan[2] >= tspan[1]`
with the explicit Dormand-Prince 5(4) Runge-Kutta method (DOPRI5) and adaptive steps. It is meant for small,
non-stiff problems whose state `y` is a `Float64` or an `SVector` of `Float64`, for which it does not allocate
inside the step loop.

The local error estimate of each step is compared with `abstol + reltol * max(|y|, |y_new|)` in each component,
in the root-mean-square norm (as in OrdinaryDiffEq.jl and Hairer's DOPRI5). The step size is chosen by Hairer's
initial-step heuristic and then by his PI controller (safety factor 0.9, step ratio between 0.2 and 10). The fifth-
order solution is propagated and the last stage is reused as the first stage of the next step, except after a stop.

Every time in `tstops` that lies inside `tspan` is stepped onto exactly; put the times where `f` is not smooth
or jumps (e.g. the ends of mission segments) there. Repeated times and times outside `tspan` are ignored. Each stop
and the ends of `tspan` are treated as possible jumps of `f`: a step that ends on a stop evaluates `f` just before
it (`prevfloat`), and the step that starts from it evaluates `f` just after it (`nextfloat`), so `f` is only used
on the open intervals between stops and its value exactly at a stop does not matter. The solution is saved at
`tspan[1]` and after every accepted step, so `ts` contains every stop and ends at `tspan[2]`.

The integration stops early, with a warning, if the step size underflows or after `maxiters` steps.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `f`: derivative function, `f(y, t)` returns a value of the same type as `y`.
    - `y0`: initial state, a `Float64` or an `SVector{N, Float64}`.
    - `tspan::Tuple{Float64, Float64}`: start and end of the integration.
    - `tstops`: times that the integration must step onto.
    - `reltol::Float64`, `abstol::Float64`: relative and absolute tolerances on the local error of each step.
    - `maxiters::Int`: maximum number of steps, accepted or rejected.

    **Outputs:**
    - `ts::Vector{Float64}`: times at which the solution is saved.
    - `ys::Vector`: solution at those times.
    - `success::Bool`: true if the integration reached `tspan[2]`.
"""
function ode_solve(f, y0::S, tspan::Tuple{Float64, Float64}; tstops = Float64[], reltol::Float64 = 1e-6,
                   abstol::Float64 = 1e-6, maxiters::Int = 100_000) where {S <: Union{Float64, SVector{<:Any, Float64}}}
    t0, tf = tspan
    tf >= t0 || throw(ArgumentError("ode_solve integrates forward in time only, got tspan = $tspan"))
    ts = [t0]
    ys = [y0]
    tf == t0 && return ts, ys, true

    # Times to step onto, ending with tf
    stops = sort!(unique!(Float64[t for t in tstops if t0 < t < tf]))
    push!(stops, tf)

    # Dormand-Prince 5(4) coefficients (Hairer, Nørsett & Wanner, Solving ODEs I, Table 5.2)
    a21 = 1 / 5
    a31, a32 = 3 / 40, 9 / 40
    a41, a42, a43 = 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63, a64, a65 = 9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656
    b1, b3, b4, b5, b6 = 35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84 # fifth order, = row 7 of a
    c2, c3, c4, c5 = 1 / 5, 3 / 10, 4 / 5, 8 / 9
    # Fifth- minus fourth-order weights
    e1, e3, e4, e5, e6, e7 = 71 / 57600, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40

    # Step-size control (Hairer's DOPRI5 defaults)
    safe = 0.9
    facmin, facmax = 0.2, 10.0 # bounds on h_new / h
    β = 0.04
    expo = 0.2 - 0.75 * β

    errnorm(e, y, yn) = sqrt(sum(abs2, e ./ (abstol .+ reltol .* max.(abs.(y), abs.(yn)))) / length(y))

    t = t0
    y = y0
    k1 = f(y, nextfloat(t))

    # Initial step (Hairer's heuristic)
    d0 = errnorm(y, y, y)
    d1 = errnorm(k1, y, y)
    h = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1
    h = min(h, tf - t0)
    k2 = f(y + h * k1, t + h)
    d2 = errnorm(k2 - k1, y, y) / h
    dmax = max(d1, d2)
    h1 = dmax <= 1e-15 ? max(1e-6, h * 1e-3) : (0.01 / dmax)^(1 / 5)
    h = min(100 * h, h1, tf - t0)

    facold = 1e-4
    rejected = false
    istop = 1
    for _ in 1:maxiters
        # Step onto the next stop if it is within reach, stretching the step by up to 1% rather than leaving a
        # sliver before the stop (as in Hairer's DOPRI5)
        tstop = stops[istop]
        lands = t + 1.01 * h >= tstop
        if lands
            h = tstop - t
        end
        tn = lands ? tstop : t + h
        te = lands ? prevfloat(tstop) : tn # time of the stages at the end of the step
        if h <= 16 * eps(abs(t)) # step size underflow
            @warn "ode_solve: step size underflow at t = $t; the integration stopped before t = $tf"
            return ts, ys, false
        end

        k2 = f(y + h * a21 * k1, t + c2 * h)
        k3 = f(y + h * (a31 * k1 + a32 * k2), t + c3 * h)
        k4 = f(y + h * (a41 * k1 + a42 * k2 + a43 * k3), t + c4 * h)
        k5 = f(y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), t + c5 * h)
        k6 = f(y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), te)
        yn = y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
        k7 = f(yn, te)
        err = errnorm(h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7), y, yn)
        isfinite(err) || (err = Inf)

        fac11 = err^expo
        if err <= 1 # accept
            fac = clamp(fac11 / facold^β / safe, 1 / facmax, 1 / facmin)
            hn = h / fac
            facold = max(err, 1e-4)
            rejected && (hn = min(hn, h)) # no growth right after a rejection
            rejected = false
            t, y, k1 = tn, yn, k7
            push!(ts, t)
            push!(ys, y)
            if lands
                t == tf && return ts, ys, true
                istop += 1
                k1 = f(y, nextfloat(t)) # f may jump at the stop
            end
            h = hn
        else # reject and retry with a smaller step
            h = h / min(1 / facmin, fac11 / safe)
            rejected = true
        end
    end
    @warn "ode_solve: maximum number of steps ($maxiters) reached at t = $t; the integration stopped before t = $tf"
    return ts, ys, false
end
