# Benchmark of the gas-property primitives: IdealGasThermo.jl-based
# `TASOPT.engine` vs the legacy table code vendored in test/legacy_thermo/.
# Not part of runtests.jl. Run with
#   julia --project -e 'include("test/benchmark_thermo.jl")'
# (BenchmarkTools must be loadable, e.g. from the test target or the
# global environment.)
using TASOPT, BenchmarkTools, Printf
isdefined(@__MODULE__, :LegacyThermo) ||
    include(joinpath(@__DIR__, "legacy_thermo", "LegacyThermo.jl"))

const E = TASOPT.engine
const L = LegacyThermo

air = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
n = length(air)
lam = (air .+ 0.025 .* E.gasfuel(24, n + 1)[1:n]) ./ 1.025
air6 = [air; 0.0]
beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

function state(M, a, T)
    s, _, h, _, cp, r = M.gassum(a, length(a), T)
    return s, h, cp, r
end

cases = [
    ("gasfun(N2, 1200 K)", M -> (() -> M.gasfun(1, 1200.0))),
    ("gassum(air, 5 species)", M -> (() -> M.gassum(air, n, 800.0))),
    ("gassumd(air, 5 species)", M -> (() -> M.gassumd(air, n, 800.0))),
    ("gas_tset(air, h(900 K))", M -> (h = M.gassum(air, n, 900.0)[3]; () -> M.gas_tset(air, n, h, 500.0))),
    ("gas_prat(air, PR=10)", M -> ((s, h, cp, r) = state(M, air, 400.0); () -> M.gas_prat(air, n, 1e5, 400.0, h, s, cp, r, 10.0, 0.9))),
    ("gas_delh(burned, -4e5)", M -> ((s, h, cp, r) = state(M, lam, 1600.0); () -> M.gas_delh(lam, n, 1e6, 1600.0, h, s, cp, r, -4e5, 1 / 0.9))),
    ("gas_burn(C14H30)", M -> (g = M.gasfuel(24, 6); () -> M.gas_burn(air6, beta, g, 6, 24, 800.0, 300.0, 1600.0, 0.0))),
]

println(rpad("function", 26), lpad("legacy", 12), lpad("new", 12), lpad("speedup", 10))
for (name, mk) in cases
    f_old = mk(L)
    f_new = mk(E)
    t_old = median(@benchmark $f_old()).time
    t_new = median(@benchmark $f_new()).time
    @printf("%-26s %9.1f ns %9.1f ns %8.2fx\n", name, t_old, t_new, t_old / t_new)
end
