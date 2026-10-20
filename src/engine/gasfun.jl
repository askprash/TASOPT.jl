using Zygote

include("gasdata.jl")

"""
    gasfun(igas, t)

Computes properties of a thermally-perfect gas
with some variable specific heat cp[T].

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `igas`: index specifying the gas (see if blocks below for list)
    - `t`: temperature T in Kelvin

    **Output:**
    - `s`: entropy-complement function s[T]
    - `s_t`: ds/dT
    - `h`: complete enthalpy function h[T]
    - `h_t`: dh/dT
    - `cp`: specific heat cp[T]
    - `r`: ideal-gas constant R


The adiabatic pressure change over a process 1->2 with some polytropic efficiency epol is:

``\\ p2  = \\ p1  exp [   epol   (s2-s1)/R ] ``   compression

``\\ p2  = \\ p1  exp [ (1/epol) (s2-s1)/R ] ``   expansion
"""
function gasfun(igas, t)

    if (igas == 1)
        s, h, cp, r = gas_N2(t)
    elseif (igas == 2)
        s, h, cp, r = gas_O2(t)
    elseif (igas == 3)
        s, h, cp, r = gas_CO2(t)
    elseif (igas == 4)
        s, h, cp, r = gas_H2O(t)
    elseif (igas == 5)
        s, h, cp, r = gas_Ar(t)
    elseif (igas == 11)
        s, h, cp, r = gas_CH4(t)
    elseif (igas == 12)
        s, h, cp, r = gas_C2H6(t)
    elseif (igas == 13)
        s, h, cp, r = gas_C3H8(t)
    elseif (igas == 14)
        s, h, cp, r = gas_C4H10(t)
    elseif (igas == 18)
        s, h, cp, r = gas_C8H18(t)
    elseif (igas == 24)
        s, h, cp, r = gas_C14H30(t)
    elseif (igas == 40)
        s, h, cp, r = gas_H2(t)
    else
        error("GASFUN: undefined gas index: ", igas, " & code=0")
    end
    s_t = cp / t
    h_t = cp
    return s, s_t, h, h_t, cp, r
end # gasfun

"""
    gaschem(igas)

Returns number of C,H,O,N atoms in gas molecule, for the gases implemented in function gasfun above.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `igas`: index specifying the gas (see if blocks below for list)

    **Output:**
    - `nchon(.)`: number of C,H,O,N atoms in gas molecule

"""
function gaschem(igas)

    kc = 1
    kh = 2
    ko = 3
    kn = 4

    nchon = [0, 0, 0, 0]
    buf = Zygote.Buffer(nchon, length(nchon))
    for i = 1:length(nchon)
        buf[i] = nchon[i]
    end

    if (igas == 1)  #   N2
        buf[kn] = 2
    elseif (igas == 2)  #   O2
        buf[ko] = 2
    elseif (igas == 3)  #   CO2
        buf[kc] = 1
        buf[ko] = 2
    elseif (igas == 4)  #   H2O
        buf[kh] = 2
        buf[ko] = 1
    elseif (igas == 5)  #   Ar

    elseif (igas == 11)  #   CH4
        buf[kc] = 1
        buf[kh] = 4
    elseif (igas == 12)  #   C2H6
        buf[kc] = 2
        buf[kh] = 6
    elseif (igas == 13)  #   C3H8
        buf[kc] = 3
        buf[kh] = 8
    elseif (igas == 14)  #   C4H10
        buf[kc] = 4
        buf[kh] = 10
    elseif (igas == 18)  #   C8H18
        buf[kc] = 8
        buf[kh] = 18
    elseif (igas == 24)  #   C14H30 Jet-A surrogate used by the legacy sizing model
        buf[kc] = 14
        buf[kh] = 30
    elseif (igas == 40)  #   H2
        buf[kc] = 0
        buf[kh] = 2
    else
        error("GASFUN: undefined gas index: ", igas, " & code=0")
    end
    nchon = copy(buf)
    return nchon
end # gaschem

function gas_N2(t1, t, tl, cp, cpt, h, s)
    r = 296.94
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_N2

function gas_O2(t1, t, tl, cp, cpt, h, s)
    r = 259.82
    hform = 0.000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_O2

function gas_Ar(t1, t, tl, cp, cpt, h, s)
    r = 208.00
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_Ar

function gas_CO2(t1, t, tl, cp, cpt, h, s)
    r = 188.96
    hform = -0.89430E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_CO2

function gas_H2O(t1, t, tl, cp, cpt, h, s)
    r = 461.91
    hform = -0.13433E+08
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_H2O

function gas_CH4(t1, t, tl, cp, cpt, h, s)
    r = 519.65
    hform = -0.46750E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_CH4

function gas_C2H6(t1, t, tl, cp, cpt, h, s)
    r = 277.15
    hform = -0.27930E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C2H6

function gas_C3H8(t1, t, tl, cp, cpt, h, s)
    r = 188.50
    hform = -0.23570E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C3H8

function gas_C4H10(t1, t, tl, cp, cpt, h, s)
    r = 143.30
    hform = -0.21640E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C4H10

function gas_C8H18(t1, t, tl, cp, cpt, h, s)
    r = 72.900
    hform = -0.18280E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C8H18

function gas_C14H30(t1, t, tl, cp, cpt, h, s)
    r = 167.00
    hform = -0.16700E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
    # gas_C14H30 
end

function gas_H2(t1, t, tl, cp, cpt, h, s)
    r = 4124.9
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_H2

"""
    findsegment(x::T, xarr::Vector{T})

Uses bisection to find the right interval of the array `xarr` where x lies.
Returns im and io s.t. `xarr[im] < x < xarr[io]`.

Additionally returns the interval `dx = xarr[io] - xarr[im]`
"""
function findsegment(x::Float64, xarr::AbstractArray)

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end

    io::Int = length(xarr)

    if x ≤ xarr[1]
        im = 1
        io = 2
    elseif x ≥ xarr[end]
        im = io - 1
    else
        im = searchsortedlast(xarr, x)
        io = im + 1
    end

    dx = xarr[io] - xarr[im]
    t = (x - xarr[im]) / dx

    return io, im, dx, t
end

"""
Helper function that returns the thermo properties by table interpolation
"""
function get_thermo(
    r::T,
    hform::T,
    t1::T,
    t::Vector{T},
    tl::Vector{T},
    cp::Vector{T},
    cpt::Vector{T},
    h::Vector{T},
    s::Vector{T},
) where {T<:AbstractFloat}

    i, im, dt, f = findsegment(t1, t)

    dtl = tl[i] - tl[im]
    fl = (log(t1) - tl[im]) / dtl
    A = 1.0 - f
    B = 1.0 - fl
    #- - - - - - - - - - -
    s1 =
        B * s[im] +
        fl * s[i] +
        fl * B * (B * (dtl * cp[im] - s[i] + s[im]) - fl * (dtl * cp[i] - s[i] + s[im]))

    h1 =
        hform +
        A * h[im] +
        f * h[i] +
        f * A * (A * (dt * cp[im] - h[i] + h[im]) - f * (dt * cp[i] - h[i] + h[im]))

    cp1 =
        A * cp[im] +
        f * cp[i] +
        f * A * (A * (dt * cpt[im] - cp[i] + cp[im]) - f * (dt * cpt[i] - cp[i] + cp[im]))

    r1 = r
    #- - - - - - - - - - -
    return s1, h1, cp1, r1
end  # function get_thermo

const _IGT_N2 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("N2"))
const _IGT_O2 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("O2"))
const _IGT_Ar = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("Ar"))
const _IGT_CO2 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("CO2"))
const _IGT_H2O = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("H2O"))
const _IGT_CH4 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("CH4"))
const _IGT_C2H6 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("C2H6"))
const _IGT_C3H8 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("C3H8"))
const _IGT_C4H10 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("C4H10,n-butan"))
const _IGT_C8H18 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("C8H18,n-octan"))
const _IGT_JETA = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("Jet-A(g)"))
const _IGT_H2 = IdealGasThermo.FrozenGas(IdealGasThermo.species_in_spdict("H2"))

@inline function _fast_thermo(gas, t)
    props = IdealGasThermo.props(gas, t)
    return props.s0, props.h, props.cp, IdealGasThermo.R(gas)
end

@inline function _cp_T(gas, t)
    a = IdealGasThermo.coeffs(gas, t)
    tinv = inv(t)
    tinv2 = tinv * tinv
    return IdealGasThermo.Runiv * (
        -2 * a[1] * tinv2 * tinv - a[2] * tinv2 +
        a[4] + t * (2 * a[5] + t * (3 * a[6] + t * 4 * a[7]))
    )
end

@inline function _fast_thermod(gas, t)
    props = IdealGasThermo.props(gas, t)
    cp = props.cp
    return props.s0, cp / t, props.h, cp, cp, _cp_T(gas, t), IdealGasThermo.R(gas)
end

# Lower temperature limit of the NASA-9 fits in IdealGasThermo.
const _IGT_TMIN = 200.0

# H2 below _IGT_TMIN (LH2 coolant, cold H2 fuel): use the tabulated data in
# gasdata.jl (20.369-200 K), shifted so that h and s equal the NASA-9 values
# at _IGT_TMIN. h and s are then continuous at the switch; cp jumps by 0.6%.
const _H2_TAB_TMIN = get_thermo(4124.9, 0.0, _IGT_TMIN, H2.t, H2.tl, H2.cp, H2.cpt, H2.h, H2.s)
const _H2_IGT_TMIN = IdealGasThermo.props(_IGT_H2, _IGT_TMIN)

@inline function _h2_table(t)
    s, h, cp, _ = get_thermo(4124.9, 0.0, t, H2.t, H2.tl, H2.cp, H2.cpt, H2.h, H2.s)
    return s - _H2_TAB_TMIN[1] + _H2_IGT_TMIN.s0, h - _H2_TAB_TMIN[2] + _H2_IGT_TMIN.h, cp
end

function _h2_cryo_thermo(t)
    s, h, cp = _h2_table(t)
    return s, h, cp, IdealGasThermo.R(_IGT_H2)
end

function _h2_cryo_thermod(t)
    s, h, cp = _h2_table(t)
    dt = 0.01
    cp_t = (_h2_table(t + dt)[3] - _h2_table(t - dt)[3]) / (2 * dt)
    return s, cp / t, h, cp, cp, cp_t, IdealGasThermo.R(_IGT_H2)
end

function gasfund(igas, t)
    if (igas == 1)
        return _fast_thermod(_IGT_N2, t)
    elseif (igas == 2)
        return _fast_thermod(_IGT_O2, t)
    elseif (igas == 3)
        return _fast_thermod(_IGT_CO2, t)
    elseif (igas == 4)
        return _fast_thermod(_IGT_H2O, t)
    elseif (igas == 5)
        return _fast_thermod(_IGT_Ar, t)
    elseif (igas == 11)
        return _fast_thermod(_IGT_CH4, t)
    elseif (igas == 12)
        return _fast_thermod(_IGT_C2H6, t)
    elseif (igas == 13)
        return _fast_thermod(_IGT_C3H8, t)
    elseif (igas == 14)
        return _fast_thermod(_IGT_C4H10, t)
    elseif (igas == 18)
        return _fast_thermod(_IGT_C8H18, t)
    elseif (igas == 24)
        return _fast_thermod(_IGT_JETA, t)
    elseif (igas == 40)
        return t < _IGT_TMIN ? _h2_cryo_thermod(t) : _fast_thermod(_IGT_H2, t)
    else
        error("GASFUN: undefined gas index: ", igas, " & code=0")
    end
end

gas_N2(t1) = _fast_thermo(_IGT_N2, t1)
gas_O2(t1) = _fast_thermo(_IGT_O2, t1)
gas_Ar(t1) = _fast_thermo(_IGT_Ar, t1)
gas_CO2(t1) = _fast_thermo(_IGT_CO2, t1)
gas_H2O(t1) = _fast_thermo(_IGT_H2O, t1)
gas_CH4(t1) = _fast_thermo(_IGT_CH4, t1)
gas_C2H6(t1) = _fast_thermo(_IGT_C2H6, t1)
gas_C3H8(t1) = _fast_thermo(_IGT_C3H8, t1)
gas_C4H10(t1) = _fast_thermo(_IGT_C4H10, t1)
gas_C8H18(t1) = _fast_thermo(_IGT_C8H18, t1)
gas_C14H30(t1) = _fast_thermo(_IGT_JETA, t1)
gas_H2(t1) = t1 < _IGT_TMIN ? _h2_cryo_thermo(t1) : _fast_thermo(_IGT_H2, t1)
