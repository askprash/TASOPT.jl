using DocStringExtensions
abstract type AbstractLayout end

"""
$TYPEDEF

Fuselage Layout Structure:
Contains dimensions, heights, etc. to design a fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct FuselageLayout <: AbstractLayout
    """Cross section definition"""
    cross_section::AbstractCrossSection = SingleBubble()
    """Thickness of webs """
    thickness_webs::Float64 = 0 #nfwebs
    """X position of nose [m]"""
    x_nose::Float64  = 0# = ac.parg[igxnose] #xnose
    """X position of pressure shell forward [m]"""
    x_pressure_shell_fwd::Float64 = 0# = ac.parg[igxshell1] #xshell1
    """X position of pressure shell aft [m]"""
    x_pressure_shell_aft::Float64 = 0# = ac.parg[igxshell2] #xshell2
    """X position of cylinder start [m]"""
    x_start_cylinder::Float64 = 0# = ac.parg[igxblend1] #xblend1
    """X position of cylinder end [m]"""
    x_end_cylinder::Float64 = 0# = ac.parg[igxblend2] #xblend2
    """X position of fuselage cone end [m]"""
    x_cone_end::Float64 = 0# = ac.parg[igxend] #xconeend
    """X position of fuselage end [m]"""
    x_end::Float64 = 0# = ac.parg[igxend] #xend
    """Tailcone taper (lambdac) [m]"""
    taper_tailcone::Float64 = 0# lambdac
    """Floor depth (depth of floor beams) [m]"""
    floor_depth::Float64 = 0
    """Nose Radius [m]"""
    nose_radius::Float64 = 0
    """Tail Radius [m]"""
    tail_radius::Float64 = 0
    """Taper fuselage to Point (0) or Edge (1)"""
    taper_fuse::Int64 = 0 # 0 = point ; 1 = edge
end

# Helper function to be able to simplify 
function Base.getproperty(layout::FuselageLayout, sym::Symbol)
    cross_section = getfield(layout, :cross_section)
    
    if sym === :l_nose
        return getfield(layout, :x_pressure_shell_fwd) - getfield(layout, :x_nose)
    elseif sym === :l_shell
        return getfield(layout, :x_pressure_shell_aft) - 
               getfield(layout, :x_pressure_shell_fwd)
    elseif sym === :l_floor
        return getproperty(layout, :l_shell) + 2.0*getproperty(layout, :radius)
    
    elseif sym ∈ (:radius, :n_webs, :bubble_lower_downward_shift, :bubble_center_y_offset)
        return getproperty(cross_section, sym)
    else
        return getfield(layout, sym)
    end
end

"""
$TYPEDEF

Wing Layout Structure:
Contains dimensions, heights, etc. to design a Wing

$TYPEDFIELDS
"""
@kwdef mutable struct WingLayout
    """Aspect Ratio [m]"""
    AR::Float64 = 0
    """Sweep [degrees]"""
    sweep::Float64 = 0
    """Wing Span [m]"""
    b::Float64 = 0
    """Span of inner wing (break/"snag") [m]"""
    b_inner::Float64 = 0
    """Max Wing Span [m]"""
    b_max::Float64 = 0
    """Outer or "tip" taper ratio of chord"""
    λt::Float64 = 0
    """Inner or break/"snag" taper ratio of chord"""
    λs::Float64 = 0
    """Span fraction of inner wing break ("snag")"""
    ηs::Float64 = 0
    """Wing center box width [m]"""
    box_width::Float64 = 0
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0
end