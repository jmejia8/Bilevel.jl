module LowerLevelSABO
import Statistics: var

import Metaheuristics
################################################################################
#                       L o w e r      L e v e l
################################################################################


function Metaheuristics.stop_criteria!(
        status,
        parameters::Metaheuristics.ECA,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    status.stop = var(map(s -> s.f, status.population)) < 1e-5

end

function initialize_eca!(
        parameters::Metaheuristics.ECA,
        problem,
        information,
        options,
        args...;
        bilevel_population = nothing,
        kargs...
    )
#=
function initialize_eca!(
    problem,
    engine,
    parameters,
    status,
    information,
    options;
    bilevel_population = nothing,
)
=#
    if isnothing(bilevel_population)
        return Metaheuristics.gen_initial_state(problem,parameters,information,options)
    end

    # parameters.N -= length(bilevel_population)


    population = []
    for i = 1:length(bilevel_population)
        sol = bilevel_population[i]
        push!(
            population,
            Metaheuristics.generateChild(sol.y, problem.f(sol.y)),
        )
    end

    parameters.N = length(status.population)

    best_sol = Metaheuristics.get_best(population)
    status = State(best_sol, population)
    status.f_calls = parameters.N
    status.iteration = 1

    return status
end

function gent_optimal_SABO(
    f,
    bounds,
    local_population = nothing;
    f_calls_limit = isnothing(local_population) ? 1000 * size(bounds, 2) :
                    100 * size(bounds, 2),
)
    D = size(bounds, 2)
    K = isnothing(local_population) ? 7 : 3
    η_max = 1.2 #isnothing(local_population) ? 1.2 : 1.2

    eca = Metaheuristics.ECA(K = K, N = K * size(bounds, 2), η_max = η_max)
    if !isnothing(local_population)

        population = []
        for i = 1:length(local_population)
            sol = local_population[i]
            push!(
                  population,
                  Metaheuristics.generateChild(sol.y, f(sol.y)),
                 )
        end

        eca.parameters.N = length(population)
        eca.status = Metaheuristics.State(Metaheuristics.get_best(population), population)

    end

    res = Metaheuristics.optimize(f, bounds, eca)


    return res
end


end
