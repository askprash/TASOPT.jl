"""
    LegacyThermo

Frozen reference copy of TASOPT's original table-based gas thermodynamics,
used only by the test suite to check the IdealGasThermo.jl-based
implementation in `src/engine` against the code it replaced.

`gasfun.jl`, `gascalc.jl` and `gasdata.jl` in this folder are verbatim copies
of `src/engine/` at commit b9be1deb (upstream `main` before the port). Do not
edit them; they are the reference, not production code.
"""
module LegacyThermo
using Zygote
include("gasfun.jl")
include("gascalc.jl")
end
