#####################################################
#
#         STRUCTURES FOR THE SOLUTIONS
#
#####################################################

mutable struct xf_indiv # Single Objective
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
end

mutable struct xfg_indiv # Single Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
end

mutable struct xfgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end

mutable struct xFgh_indiv # multi Objective Constraied
    x::Vector{Float64}
    y::Vector{Float64}
    F::Vector{Float64}
    f::Vector{Float64}
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end


#####################################################
#
#         STRUCTURES FOR THE OPTIONS
#
#####################################################

struct Options
    # upper level parameters
    x_tol::Real
    F_tol::Real
    G_tol::Real
    H_tol::Real
    F_calls_limit::Int
    G_calls_limit::Int
    H_calls_limit::Int
    
    # lower level parameters
    y_tol::Real
    f_tol::Real
    g_tol::Real
    h_tol::Real
    f_calls_limit::Int
    g_calls_limit::Int
    h_calls_limit::Int
    
    iterations::Int
    store_convergence::Bool
    show_results::Bool
    debug::Bool
end

function Options(;
    # upper level parameters
    x_tol::Real = 0.0,
    F_tol::Real = 0.0,
    G_tol::Real = 0.0,
    H_tol::Real = 0.0,
    F_calls_limit::Int = 0,
    G_calls_limit::Int = 0,
    H_calls_limit::Int = 0,
    
    # lower level parameters
    y_tol::Real = 0.0,
    f_tol::Real = 0.0,
    g_tol::Real = 0.0,
    h_tol::Real = 0.0,
    f_calls_limit::Int = 0,
    g_calls_limit::Int = 0,
    h_calls_limit::Int = 0,
    
    iterations::Int = 1000,
    store_convergence::Bool = false,
    show_results::Bool = true,
    debug::Bool = false)

    
    Options(
        # upper level parameters
        promote(x_tol, F_tol, G_tol, H_tol)...,
        promote(F_calls_limit, G_calls_limit, H_calls_limit)...,
        
        # lower level parameters
        promote(y_tol, f_tol, g_tol, h_tol)...,
        promote(f_calls_limit, g_calls_limit, h_calls_limit)...,
        
        Int(iterations),

        # Results options
        promote(store_convergence,show_results, debug)...
    )
    
end

#####################################################
#
#         STRUCTURES FOR THE RESULTS
#
#####################################################

mutable struct Results
    # upper level parameters
    Δx::Real
    ΔF::Real
    ΔG::Real
    ΔH::Real
    F_calls::Int
    G_calls::Int
    H_calls::Int
    
    # lower level parameters
    Δy::Real
    Δf::Real
    Δg::Real
    Δh::Real
    f_calls::Int
    g_calls::Int
    h_calls::Int
    
    iterations::Int
    best_sol::DataType
    # convergence::
end

#####################################################
#
#         STRUCTURES FOR THE ITERATION STATE
#
#####################################################

mutable struct State{T<:Int}
    best_sol::DataType

    # upper level parameters
    F_calls::T
    G_calls::T
    H_calls::T

    # upper level parameters
    f_calls::T
    g_calls::T
    h_calls::T

    iteration::T

end