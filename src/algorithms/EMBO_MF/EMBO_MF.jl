mutable struct MBO
    K_ul::Int
    K_ll::Int
    N::Int
    η_max_ul::Float64
    η_max_ll::Float64
    λ_ul::Array{Array{Float64}}
    λ_ll::Array{Array{Float64}}
    w::Vector{Float64}
    z::Vector{Float64}
    B_ul::Array{Array{Int}}
    B_ll::Array{Array{Int}}
    T_ul::Int
    T_ll::Int
    Ψ::Function
    G_te::Vector{Float64}
end


function update_neighbors_ul!(parameters)
    N = length(parameters.λ_ul)
    d = zeros(N, N)

    empty!(parameters.B_ul)

    for i in 1:N
        for j in (i+1):N
            d[i, j] = norm(parameters.λ_ul[i] - parameters.λ_ul[j])
            d[j, i] = d[i, j]
        end

        I = sortperm(d[i,:])
        push!(parameters.B_ul, I[2:(parameters.T_ul+1)])
    end

    parameters.B_ul
end

function initialize_MBO!(problem, engine, parameters, status, information, options)
  
    status.final_time = time()
end



function update_state_MBO!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    t,
)
 

end


function is_better_MBO(
    New::xFgh_indiv,
    Old::xFgh_indiv;
    ε = 0.0,
)

    old_vio = Metaheuristics.violationsSum(Old.G, Old.H, ε = ε) # upper level
    old_vio += Metaheuristics.violationsSum(Old.g, Old.h, ε = ε) # lower level

    new_vio = Metaheuristics.violationsSum(New.G, New.H, ε = ε) # upper level
    new_vio += Metaheuristics.violationsSum(New.g, New.h, ε = ε) # lower level

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    for i in 1:length(Old.F)
        if Old.F[i] < New.F[i]
            return false
        end
    end

    return true
end


function stop_criteria_MBO(status::State, information::Information, options::Options)
    Bilevel.ull_call_limit_stop_check(status, information, options) ||
    Bilevel.iteration_stop_check(status, information, options)
end

function final_stage_MBO!(status, information, options)
    status.final_time = time()

    # compute Pareto front if it is a multiobjective problem
    if typeof(status.population[1].F) <: Array && isempty(status.best_sol)
        options.debug && @info "Computing Pareto front..."
        status.best_sol = Metaheuristics.get_pareto_front(status.population, is_better_MBO)
    end
end




function MBO(;
    K_ul = 3,
    K_ll = 3,
    N = 0,
    η_max_ul = 1.5,
    η_max_ll = 1.5,
    λ_ul = zeros(0),
    λ_ll = zeros(0),
    w = zeros(0),
    z = zeros(0),
    B_ul = [Int[]],
    B_ll = [Int[]],
    T_ul = 10,
    T_ll = 10,
    Ψ = identity,
    G_te = zeros(0),
    F_calls_limit = 1000,
    f_calls_limit = Inf,
    iterations = 1000,
    options = nothing,
    information = Information(),
)

    if isnothing(options)
        options = Bilevel.Options(
            F_calls_limit = F_calls_limit,
            f_calls_limit = Inf,
            debug = false,
            iterations = iterations,
            store_convergence = false,
        )
    end



    parameters = MBO(K_ul,
                    K_ll,
                    N,
                    η_max_ul,
                    η_max_ll,
                    λ_ul,
                    λ_ll,
                    w,
                    z,
                    B_ul,
                    B_ll,
                    T_ul,
                    T_ll,
                    Ψ,
                    G_te
                )

    algorithm = Bilevel.Algorithm(
        parameters;
        initialize! = initialize_MBO!,
        update_state! = update_state_MBO!,
        lower_level_optimizer = lower_level_optimizer_MBO,
        is_better = is_better_MBO,
        stop_criteria = stop_criteria_MBO,
        final_stage! = final_stage_MBO!,
        options = options,
        information = information,
    )



    algorithm
end


export MBO, optimize
