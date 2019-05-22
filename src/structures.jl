#####################################################
#
#         STRUCTURES FOR THE SOLUTIONS
#
#####################################################

mutable struct xf_indiv <: AbstractSolution # Single Objective
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

mutable struct Options
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
    ll_iterations::Int
    store_convergence::Bool
    show_results::Bool
    debug::Bool
    search_type::Symbol
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
    ll_iterations::Int = 1000,
    store_convergence::Bool = false,
    show_results::Bool = true,
    debug::Bool = false,
    search_type::Symbol=:minimize)

    
    Options(
        # upper level parameters
        promote(x_tol, F_tol, G_tol, H_tol)...,
        promote(F_calls_limit, G_calls_limit, H_calls_limit)...,
        
        # lower level parameters
        promote(y_tol, f_tol, g_tol, h_tol)...,
        promote(f_calls_limit, g_calls_limit, h_calls_limit)...,
        
        promote(iterations,ll_iterations)...,

        # Results options
        promote(store_convergence,show_results, debug)...,
        Symbol(search_type)
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
    population::Array

    # upper level parameters
    F_calls::T
    G_calls::T
    H_calls::T

    # upper level parameters
    f_calls::T
    g_calls::T
    h_calls::T

    iteration::T
    success_rate::Real
    convergence::Array{State}

end

function State(
        best_sol,
        population::Array;

        # upper level parameters
        F_calls::Int = 0,
        G_calls::Int = 0,
        H_calls::Int = 0,

        # upper level parameters
        f_calls::Int = 0,
        g_calls::Int = 0,
        h_calls::Int = 0,
        
        iteration::Int = 0,

        success_rate::Real = 0,
        convergence::Array{State} = State[],
    )

    State(#
        best_sol,
        Array(population),
        
        # upper level parameters
        promote(
            F_calls,
            G_calls,
            H_calls,
            f_calls,
            g_calls,
            h_calls,
            iteration)...,
            Real(success_rate),
            State[])
    
end

struct Problem
    F::Function
    f::Function
    bounds_ul::Array
    bounds_ll::Array
    G::Function
    g::Function
    type::Symbol
end

function Problem(F::Function,
                f::Function,
                bounds_ul::Array,
                bounds_ll::Array;
                G::Function = _1_(x,y) = 0,
                g::Function = _2_(x,y) = 0)

    
    type::Symbol = :constrained

    if nameof(G) == :_1_ && nameof(G) == :_2_
        type = :unconstrained
    elseif nameof(G) == :_1_
        type = :constrained_ll
    elseif nameof(G) == :_2_
        type = :constrained_ul
    end
    
    Problem(F, f, bounds_ul, bounds_ll, G, g, type)
end

#####################################################
#
#         STRUCTURES FOR THE ALGORITHMS
#
#####################################################

struct QBCA
    # QBCA Options
    k::Int
    N::Int
    η_ul_max::Real
    η_ll_max::Real
    α::Real
    β::Real
    s_min::Real

    options::Options

end

function QBCA(D_ul;
        # QBCA parameters
        k::Int = 3,
        N::Int = 2k * D_ul,
        η_ul_max::Real = 2.0,
        η_ll_max::Real = 1.0 / η_ul_max,
        s_min::Real = 0.01,
        α::Real = 0.05,
        β::Real = 0.05,


        # general Options
        iterations::Int  = 500D_ul,
        # lower level options
        ll_iterations::Int = 1000,
        
        F_calls_limit::Real = 1000D_ul,
        f_calls_limit::Real = Inf,

        options::Options = Options(),

    )

   
    options.iterations = iterations
    options.ll_iterations = ll_iterations
    options.F_calls_limit = F_calls_limit
    options.f_calls_limit = f_calls_limit

    QBCA(#
        # general Options
        promote(k, N)...,
        promote(η_ul_max,η_ll_max,α, β, s_min)...,
        options
    )

end

struct Information
    F_optimum::Real
    f_optimum::Real

    x_optimum::Array{Real}
    y_optimum::Array{Real}
end

function Information(;#
    F_optimum::Real = NaN,
    f_optimum::Real = NaN,
    x_optimum::Array{Real} = Real[],
    y_optimum::Array{Real} = Real[],
    )

    Information(promote(F_optimum,f_optimum)..., x_optimum, y_optimum)

end

mutable struct Algorithm
    parameters
    initial_state::State
    update_state!::Function
    lower_level_optimizer::Function
    is_better::Function
    stop_criteria::Function
    final_stage!::Function
    information::Information
    options::Options
end

function Algorithm(   parameters,
                   initial_state::State    = State(nothing, []),
                      initialize!::Function = _1(kwargs...) = nothing;
                   update_state!::Function = _2(kwargs...) = nothing,
           lower_level_optimizer::Function = _3(kwargs...) = nothing,
                       is_better::Function = is_better, # is_better(a, b)  = true if x is better that y 
                   stop_criteria::Function = stop_check,
                    final_stage!::Function = _4(kwargs...) = nothing,
                     information::Information = Information(),
                         options::Options  = Options())
    


    Algorithm(  parameters,
                initial_state,
                update_state!,
                lower_level_optimizer,
                is_better,
                stop_criteria,
                final_stage!,
                information,
                options)

end
