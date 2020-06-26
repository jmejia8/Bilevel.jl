#####################################################
#
#         STRUCTURES FOR THE SOLUTIONS
#
#####################################################

mutable struct xf_indiv <: AbstractSolution # Single Objective
    x::Vector{Float64}
    y
    F::Float64
    f::Float64
end

mutable struct xfg_indiv <: AbstractSolution # Single Objective Constraied
    x::Vector{Float64}
    y
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
end

mutable struct xfgh_indiv <: AbstractSolution # Single Objective Constraied
    x::Vector{Float64}
    y
    F::Float64
    f::Float64
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end

mutable struct xFgh_indiv <: AbstractSolution # multi Objective Constraied
    x::Vector{Float64}
    y
    F::Vector{Float64}
    f::Vector{Float64}
    G::Vector{Float64}
    g::Vector{Float64}
    H::Vector{Float64}
    h::Vector{Float64}
end

function Individual()
    #function body
end


#####################################################
#
#         STRUCTURES FOR THE OPTIONS
#
#####################################################

mutable struct Options
    # upper level parameters
    x_tol::Float64
    F_tol::Float64
    G_tol::Float64
    H_tol::Float64
    F_calls_limit::Float64
    G_calls_limit::Float64
    H_calls_limit::Float64

    # lower level parameters
    y_tol::Float64
    f_tol::Float64
    g_tol::Float64
    h_tol::Float64
    f_calls_limit::Float64
    g_calls_limit::Float64
    h_calls_limit::Float64

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
        promote(Float64(x_tol), F_tol, G_tol, H_tol)...,
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
    Δx::Float64
    ΔF::Float64
    ΔG::Float64
    ΔH::Float64
    F_calls::Int
    G_calls::Int
    H_calls::Int

    # lower level parameters
    Δy::Float64
    Δf::Float64
    Δg::Float64
    Δh::Float64
    f_calls::Int
    g_calls::Int
    h_calls::Int

    iterations::Int
    best_sol
    # convergence::
end

mutable struct LLResult
    # lower level info
    y
    f
    Δy::Float64
    Δf::Float64
    Δg::Float64
    Δh::Float64
    f_calls::Int
    g_calls::Int
    h_calls::Int

    iterations::Int
    other
end

function LLResult(y,f;Δy = 0.0,
                    Δf = 0.0,
                    Δg = 0.0,
                    Δh = 0.0,
                    f_calls = 0,
                    g_calls = 0,
                    h_calls = 0,
                    iterations=0,
                    other=nothing)

    LLResult(y,f,promote(Δy,Δf,Δg,Δh)...,
              promote(f_calls,g_calls,h_calls,iterations)...,
              other)
end

#####################################################
#
#         STRUCTURES FOR THE ITERATION STATE
#
#####################################################

mutable struct State
    best_sol
    population::Array

    # upper level parameters
    F_calls::Int
    G_calls::Int
    H_calls::Int

    # upper level parameters
    f_calls::Int
    g_calls::Int
    h_calls::Int

    iteration::Int
    success_rate::Float64
    convergence::Array{State}
    initial_time::Float64
    final_time::Float64
    stop::Bool
    stop_msg::String

end

function State(
        best_sol,
        population;

        # upper level parameters
        F_calls = 0,
        G_calls = 0,
        H_calls = 0,

        # upper level parameters
        f_calls = 0,
        g_calls = 0,
        h_calls = 0,

        iteration= 0,

        success_rate= 0,
        convergence = State[],
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
            State[],
            0.0,
            0.0,
            false,
            ""
            )

end

struct Problem
    F::Function
    f::Function
    bounds_ul::Matrix{Float64}
    bounds_ll::Matrix{Float64}
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

    Problem(F, f, Matrix{Float64}(bounds_ul), Matrix{Float64}(bounds_ll), G, g, type)
end

#####################################################
#
#         STRUCTURES FOR THE ALGORITHMS
#
#####################################################

struct Information
    F_optimum::Float64
    f_optimum::Float64

    x_optimum::Array{Float64}
    y_optimum::Array{Float64}
end

function Information(;#
    F_optimum::Real = NaN,
    f_optimum::Real = NaN,
    x_optimum::Array{Real} = Real[],
    y_optimum::Array{Real} = Real[],
    )

    Information(promote(Float64(F_optimum),f_optimum)..., x_optimum, y_optimum)

end

mutable struct Engine
    initialize!::Function
    update_state!::Function
    lower_level_optimizer::Function
    is_better::Function
    stop_criteria::Function
    final_stage!::Function
end

function Engine(;initialize!::Function = _1(kwargs...) = nothing,
                   update_state!::Function = _2(kwargs...) = nothing,
           lower_level_optimizer::Function = _3(kwargs...) = nothing,
                       is_better::Function = _4(kwargs...) = false,
                   stop_criteria::Function = _5(kwargs...) = nothing,
                    final_stage!::Function = _6(kwargs...) = nothing)

    Engine(initialize!,update_state!,lower_level_optimizer,
           is_better,stop_criteria,final_stage!)
end

mutable struct Algorithm
    parameters
    status::State
    information::Information
    options::Options
    engine::Engine
end

function Algorithm(   parameters;
                   initial_state::State    = State(nothing, []),
                      initialize!::Function = _1(kwargs...) = nothing,
                   update_state!::Function = _2(kwargs...) = nothing,
           lower_level_optimizer::Function = _3(kwargs...) = nothing,
                       # is_better(a, b)  = true if x is better that y
                       is_better::Function = _5(kwargs...) = false,
                   stop_criteria::Function = stop_check,
                    final_stage!::Function = _4(kwargs...) = nothing,
                     information::Information = Information(),
                         options::Options  = Options())


    engine = Engine(initialize!,
                update_state!,
                lower_level_optimizer,
                is_better,
                stop_criteria,
                final_stage!)

    Algorithm(  parameters,
                initial_state,
                information,
                options,
                engine)

end


# function Algorithm(   parameters;
#                    initial_state::State    = State(nothing, []),
#                      information::Information = Information(),
#                          options::Options  = Options(),
#                          engine::Engine = Engine())
# end
