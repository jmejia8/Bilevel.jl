function lower_level_optimizer_MBO(x, problem, status, information, options, t; lower_level_parameters = nothing, n_followers)

    D_ll = size(problem.bounds_ll, 2)
    n_objectives_ll = 2#length(status.population.f)
    if isnothing(lower_level_parameters)
        options.debug &&
        @warn "No upper level parameters detected, using defeault parameters (slow)"
        lower_level_parameters = Metaheuristics.MOEAD_DE(D_ll, n_objectives_ll)
    end

    parameters = lower_level_parameters
    # empty!(parameters.population)
    # empty!(parameters.best_sol)


    Y = []

    f_calls = 0
    for i in 1:n_followers
        res = Metaheuristics.optimize( y -> problem.f(x, y, i), problem.bounds_ll, parameters)
        options.debug && display(res)
        f_calls += res.f_calls
        push!(Y, res )
    end

        
    return LLResult(Y, nothing; f_calls = f_calls)
end
