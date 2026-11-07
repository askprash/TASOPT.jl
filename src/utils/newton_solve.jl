using LinearAlgebra: norm

"""
    newton_solve(f, x0; ftol = 1e-8, xtol = 0.0, itmax = 100)

Solves the square nonlinear system `f(x) = 0` by Newton's method, starting from `x0`. It is meant for
small systems (a few unknowns) whose residual function `f(x)::Vector{Float64}` does not provide
derivatives.

The Jacobian is estimated by central finite differences, with the step `cbrt(eps()) * max(|x_j|, 1)`
in unknown `j` (the default of NLsolve.jl). Each iteration takes the Newton step with the largest
length `1, 1/2, 1/4, ...` (at most 20 halvings) that reduces `sum(abs2, f(x))`, so a residual function
that returns a large value where it cannot be evaluated pushes the step back. The iteration stops when
`maximum(abs, f(x)) <= ftol` or a step has a 2-norm `<= xtol` (the convergence tests of NLsolve.jl),
when no step length reduces the residual or the Jacobian is singular, and after `itmax` iterations.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `f`: residual function, `f(x::Vector{Float64})` returns a `Vector{Float64}` of the same length.
    - `x0::Vector{Float64}`: initial guess.
    - `ftol::Float64`: tolerance on the largest absolute residual.
    - `xtol::Float64`: tolerance on the 2-norm of a step.
    - `itmax::Int`: maximum number of iterations.

    **Outputs:**
    - `x::Vector{Float64}`: last iterate.
    - `converged::Bool`: true if the iteration stopped by the `ftol` or `xtol` test.
"""
function newton_solve(f, x0::AbstractVector{Float64}; ftol::Float64 = 1e-8, xtol::Float64 = 0.0,
                      itmax::Int = 100)
    n = length(x0)
    x = Vector{Float64}(x0)
    F = f(x)
    F2 = sum(abs2, F)
    J = zeros(n, n)
    xp = similar(x)
    for _ in 1:itmax
        maximum(abs, F) <= ftol && return x, true

        # Central-difference Jacobian
        for j in 1:n
            h = cbrt(eps()) * max(abs(x[j]), 1.0)
            copyto!(xp, x)
            xp[j] = x[j] + h
            Fp = f(xp)
            xp[j] = x[j] - h
            Fm = f(xp)
            for i in 1:n
                J[i, j] = (Fp[i] - Fm[i]) / (2h)
            end
        end
        dx = J \ F
        all(isfinite, dx) || return x, false # singular Jacobian

        # Largest step 1, 1/2, 1/4, ... that reduces the residual
        t = 1.0
        xn = x .- dx
        Fn = f(xn)
        F2n = sum(abs2, Fn)
        for _ in 1:20
            F2n < F2 && break
            t *= 0.5
            xn = x .- t .* dx
            Fn = f(xn)
            F2n = sum(abs2, Fn)
        end
        F2n < F2 || return x, false # no step length reduces the residual

        x, F, F2 = xn, Fn, F2n
        t * norm(dx) <= xtol && return x, true
    end
    return x, maximum(abs, F) <= ftol
end
