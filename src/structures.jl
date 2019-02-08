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
    F_calls_limit::Real
    G_calls_limit::Real
    H_calls_limit::Real
    
    # lower level parameters
    y_tol::Real
    f_tol::Real
    g_tol::Real
    h_tol::Real
    f_calls_limit::Real
    g_calls_limit::Real
    h_calls_limit::Real
    
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
    F_calls_limit::Real = 0,
    G_calls_limit::Real = 0,
    H_calls_limit::Real = 0,
    
    # lower level parameters
    y_tol::Real = 0.0,
    f_tol::Real = 0.0,
    g_tol::Real = 0.0,
    h_tol::Real = 0.0,
    f_calls_limit::Real = 0,
    g_calls_limit::Real = 0,
    h_calls_limit::Real = 0,
    
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
    best_sol

    # upper level parameters
    F_calls::T
    G_calls::T
    H_calls::T

    # upper level parameters
    f_calls::T
    g_calls::T
    h_calls::T

    iteration::T
    success_rate::T

end

function State(
        best_sol;

        # upper level parameters
        F_calls::Int = 0,
        G_calls::Int = 0,
        H_calls::Int = 0,

        # upper level parameters
        f_calls::Int = 0,
        g_calls::Int = 0,
        h_calls::Int = 0,
        
        iteration::Int = 0,

        success_rate::Int = 0,
    )

    State(#
        best_sol,
        
        # upper level parameters
        promote(
            F_calls,
            G_calls,
            H_calls,
            f_calls,
            g_calls,
            h_calls,
            iteration,
            success_rate)...)
    
end

#####################################################
#
#         STRUCTURES FOR THE ALGORITHMS
#
#####################################################

struct QBCA
    F::Function # upper level objective function
    f::Function # lower level objective function
    bounds_ul::Array
    bounds_ll::Array
    stop_criteria::Function
    search_type::Symbol

    # general Options
    k::Int
    N::Int
    η_ul_max::Real
    η_ll_max::Real
    s_min::Real

    # general Options
    iterations::Int
    F_calls_limit::Real
    f_calls_limit::Real
    options::Options

end

function QBCA(
        F::Function, # upper level objective function
        f::Function, # lower level objective function
        bounds_ul::Array,
        bounds_ll::Array;
        stop_criteria::Function = x -> false,
        search_type::Symbol = :minimize,

        # QBCA parameters
        k::Int = 3,
        N::Int = 2k * size(bounds_ul, 2),
        η_ul_max::Real = 2.0,
        η_ll_max::Real = 1.0 / η_ul_max,
        s_min::Real = 0.01,


        # general Options
        iterations::Int  = 500size(bounds_ul, 2),
        F_calls_limit::Real = 1000size(bounds_ul, 2),
        f_calls_limit::Real = Inf,

        options::Options = Options(iterations = iterations,F_calls_limit = F_calls_limit,f_calls_limit = f_calls_limit)

    )

    QBCA(#
        F, f,
        promote(bounds_ul,bounds_ll)...,
        stop_criteria,
        Symbol(search_type),

        # general Options
        promote(k, N)...,
        promote(η_ul_max,η_ll_max,s_min)...,

        Int(iterations),
        # general Options
        promote(F_calls_limit,f_calls_limit)...,
        options
    )

end