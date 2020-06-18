################################################################################
#                       L o w e r      L e v e l
################################################################################


function stop_check_eca(status, information, options)
    var(map(s -> s.f, status.population)) < 1e-5 ||
        Metaheuristics.stop_check(status, information, options)

end

function initialize_eca!(
    problem,
    engine,
    parameters,
    status,
    information,
    options;
    bilevel_population = nothing,
)
    if isnothing(bilevel_population)
        Metaheuristics.initialize_eca!(
            problem,
            engine,
            parameters,
            status,
            information,
            options,
        )
        return
    end

    parameters.N -= length(bilevel_population)
    Metaheuristics.initialize_eca!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
    )


    for i = 1:length(bilevel_population)
        sol = bilevel_population[i]
        push!(
            status.population,
            Metaheuristics.generateChild(sol.y, problem.f(sol.y)),
        )
    end

    parameters.N = length(status.population)

    status.f_calls = parameters.N
    status.iteration = 0
    status.best_sol = Metaheuristics.getBest(status.population, :minimize)
    status.stop = engine.stop_criteria(status, information, options)
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
    eca.engine.stop_criteria = stop_check_eca
    if !isnothing(local_population)
        eca.options.f_calls_limit = f_calls_limit
        eca.engine.initialize! =
            (a1, a2, a3, a4, a5, a6) -> initialize_eca!(
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                bilevel_population = local_population,
            )
    end

    res = Metaheuristics.optimize(f, bounds, eca)


    return res
end

function lower_level_optimizer_SABO(x, problem, status, information, options, t)
    D = size(problem.bounds_ll, 2)

    if length(status.population) >= 3D
        distances = map(sol -> norm(sol.x - x), status.population)
        I = sortperm(distances)
        local_population = status.population[I[1:D]]
        res = gent_optimal_SABO(
            z -> problem.f(x, z),
            problem.bounds_ll,
            local_population,
        )
    else
        res = gent_optimal_SABO(z -> problem.f(x, z), problem.bounds_ll)
    end
    y = res.best_sol.x
    fy = res.best_sol.f
    f_calls = res.f_calls


    return LLResult(y, fy; f_calls = f_calls)


end


################################################################################
#                    I n i t i l i a z a t i o n
################################################################################



function initialize_SABO!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
)
    p = problem

    a = p.bounds_ul[1, :]
    b = p.bounds_ul[2, :]

    X = a' .+ (b - a)' .* rand(parameters.N, size(p.bounds_ul, 2))

    empty!(status.population)
    best = nothing
    for i = 1:parameters.N
        x = X[i, :]
        ll_result =
            engine.lower_level_optimizer(x, p, status, information, options, 0)

        y = ll_result.y

        child = generateChild(x, y, p.F(x, y), ll_result.f)
        push!(status.population, child)

        status.f_calls += ll_result.f_calls
        status.F_calls += 1

        if best == nothing || engine.is_better(child, best)
            best = child
        end

        if options.debug
            display(child)
            println("------------------------------")
        end
    end


    status.best_sol = best
end


################################################################################
#                    I T E R A T I O N
################################################################################

function update_state_SABO!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    t_main_loop,
)
    N = length(status.population)
    I = randperm(N)
    c = zeros(size(problem.bounds_ul, 2))
    a = problem.bounds_ul[1, :]
    b = problem.bounds_ul[2, :]

    X = map(sol -> sol.x', status.population)
    y = map(sol -> sol.F, status.population)
    X = vcat(X...)
    X = (X .- a') ./ (b - a)'
    method = KernelInterpolation(y, X, λ = parameters.λ)
    train!(method)
    F̂ = approximate(method)

    # -------------------------------------------------------------
    # -------------------------------------------------------------
    # -------------------------------------------------------------
    x_initial = (status.best_sol.x - a) ./ (b - a)
    optimizer = Optim.Fminbox(Optim.BFGS())
    res = Optim.optimize(
        F̂,
        zeros(length(a)),
        ones(length(a)),
        x_initial,
        optimizer,
        Optim.Options(outer_iterations = 1),
    )
    p = a .+ (b - a) .* res.minimizer



    ll_result = engine.lower_level_optimizer(
        p,
        problem,
        status,
        information,
        options,
        t_main_loop,
    )
    status.f_calls += ll_result.f_calls
    q = ll_result.y
    FF = problem.F(p, q)
    if FF < status.best_sol.F
        status.best_sol.F = FF
        status.best_sol.f = ll_result.f
        status.best_sol.x = p
        status.best_sol.y = q

    end

    BCAOperators.update_state!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        t_main_loop,
    )

    if t_main_loop % 5 != 0 || status.stop
        return
    end
    best_y = status.best_sol.y
    sol = status.best_sol
    res = gent_optimal_SABO(z -> problem.f(sol.x, z), problem.bounds_ll, [sol])
    y = res.best_sol.x
    f_calls = res.f_calls

    FF = sol.F

    fy = res.best_sol.f
    sol.y = y
    sol.f = fy
    sol.F = problem.F(sol.x, sol.y)

    status.f_calls += f_calls
    status.F_calls += 1

    options.debug && @info "Re-evaluating solutions..."


    for sol in status.population
        res = gent_optimal_SABO(
            z -> problem.f(sol.x, z),
            problem.bounds_ll,
            [sol],
        )

        ff, FF = sol.f, sol.F

        sol.y = res.best_sol.x
        sol.f = res.best_sol.f
        FF = sol.F
        sol.F = problem.F(sol.x, sol.y)



        status.f_calls += res.f_calls
        status.F_calls += 1


        if status.best_sol.F > sol.F
            status.best_sol = sol
        end

        if FF == sol.F && ff == sol.f
            break
        end

    end
end



is_better_SABO(solution_1, solution_2) = solution_1.F < solution_2.F


function stop_criteria_SABO(status, information, options)
    if stop_check(status, information, options)
        return true
    end

    cond = var(map(s -> s.F, status.population)) < 1e-8
    options.debug && cond && @info("Stopped since ul_varF was met.")
    status.stop = cond
    if status.stop
        status.stop_msg = "ul_diversity"
    else
        status.stop_msg = ""
    end

    cond
end

function final_stage_SABO!(status, information, options) end




mutable struct SABO
    N::Int64
    K::Int64
    η_max::Float64
    λ::Float64
end

function SABO(
    D::Int;
    K = 3,
    N = K * D,
    η_max = 1.2,
    λ = 1e-5,
    F_calls_limit = 350D,
    f_calls_limit = Inf,
    iterations = 1 + round(F_calls_limit / N),
    options = nothing,
    information = Information(),
)

    if isnothing(options)
        options = Options(
            F_calls_limit = F_calls_limit,
            f_calls_limit = Inf,
            debug = false,
            F_tol = 1e-2,
            f_tol = 1e-3,
            store_convergence = true,
        )
    end

    # if isnothing(information)
    #     information = Information(F_optimum = 0.0, f_optimum = 0.0)
    # end

    sabo = SABO(N, K, η_max, λ)

    algorithm = Algorithm(
        sabo;
        initialize! = initialize_SABO!,
        update_state! = update_state_SABO!,
        lower_level_optimizer = lower_level_optimizer_SABO,
        is_better = is_better_SABO,
        stop_criteria = stop_criteria_SABO,
        final_stage! = final_stage_SABO!,
        options = options,
        information = information,
    )



    algorithm

end

optimize(F, f, bounds_ul, bounds_ll, method::SABO) =
    optimize(F, f, bounds_ul, bounds_ll, method)


# function BCA3(D; K = 3, N = K*D, η_max = 1.2, F_calls_limit = 350D, store_convergence=true)
#     options = Options(F_calls_limit  = F_calls_limit,
#                         f_calls_limit= Inf,
#                         debug=false,
#                         F_tol = 1e-2,
#                         f_tol = 1e-3,
#                         store_convergence=true)
#
#     information = Information(F_optimum = 0.0, f_optimum = 0.0)
#
#
#     method = SABOModule.SABO(D_ul; K = K,
#                                      N = N,
#                                      η_max = η_max,
#                                      options = options,
#                                      information = information)
# end
