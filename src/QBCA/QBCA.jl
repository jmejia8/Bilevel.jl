
struct QBCA
    # QBCA Options
    k::Int
    N::Int
    η_ul_max::Float64
    η_ll_max::Float64
    α::Float64
    β::Float64
    s_min::Float64
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
        options = Options(),
        information = Information(),


    )


    qbca = QBCA(#
        # general Options
        promote(k, N)...,
        promote(η_ul_max,η_ll_max,α, β, s_min)...,
    )

    algorithm = Algorithm(
        qbca;
        initialize! = initialize_QBCA!,
        update_state! = update_state_QBCA!,
        lower_level_optimizer = lower_level_optimizer_QBCA,
        is_better = is_better_SABO,
        stop_criteria = stop_criteria_SABO,
        final_stage! = final_stage_SABO!,
        options = options,
        information = information,
    )

end


function nearest(P, x; tol = 1e-5)
    x_nearest = P[1].x
    y = P[1].y
    d = Inf

    for sol in P
        n = norm(x - sol.x)

        n >= d && (continue)

        x_nearest = sol.x
        y = sol.y
        d = n

        d <= tol && (break)

    end

    y, d
end

################################################################################
################################################################################
################################################################################

function initialize_QBCA!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
)

    bounds_ul = problem.bounds_ul
    bounds_ll = problem.bounds_ll
    search_type = options.search_type


    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2)
    nevals_ul = 0


    # general parameters
    D = D_ul + D_ll
    # N = D < 5 ? 5*D : k*D

    # tolerance for y proximity
    tol = 1e-16

    # initialize population
    status.population, f_calls = init_population(problem.F, problem.f, parameters.N, bounds_ul, bounds_ll)

    status.F_calls = parameters.N
    status.f_calls += f_calls

    # current generation
    status.iteration = 0

    # best solution
    status.best_sol = getBest(status.population)



end

function update_state_QBCA!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    t_main_loop,
)

    k = parameters.k
    N = parameters.N
    η_ul_max = parameters.η_ul_max
    η_ll_max = parameters.η_ll_max
    α = parameters.α
    β = parameters.β

    I_ul = randperm(N)
    I_ll = randperm(N)

    status.success_rate = 0.0

    for i = 1:N

        # current
        x = status.population[i].x
        y = status.population[i].y

        # generate U masses
        U = getU(status.population, k, I_ul, i, N)
        V = getU(status.population, k, I_ll, i, N)

        # generate center of mass
        c_ul, c_ll, u_worst, v_worst = center(U, V, α, β, options.search_type)

        # stepsize
        η_ul = η_ul_max * rand()

        # u: worst element in U
        u = U[u_worst].x
        v = V[v_worst].y

        # current-to-center
        p = x + η_ul * (c_ul - u)
        p = correct(p, problem.bounds_ul)

        # ------------------------------------------- -- ---
        ll = engine.lower_level_optimizer(p, problem, status, information, options, t_main_loop)
        q = ll.y
        status.f_calls += ll.f_calls
        # ------------------------------------------- -- ---

        # q = correct(q, bounds_ll)

        sol = generateChild(p, q, F(p, q), fpq)
        status.F_calls += 1

        if sol ≺ status.population[i]
            status.population[getWorstInd(status.population)] = sol
            status.success_rate += 1.0 / N

            if sol ≺ status.best_sol
                status.best_sol = sol

                options.store_convergence &&
                    push!(convergence, deepcopy(status))

                # check stop condition
                stop = stop_check(status, information, options)
                stop && break
            end
        end

        stop = ull_call_limit_stop_check(status, information, options)
        stop && break

    end

    status.population = deepcopy(status.population)

    status.iteration += 1

    stop = stop || stop_check(status, information, options)



end


function lower_level_optimizer_QBCA(p, problem, status, information, options, t)
    tol = 1e-16
    y0, d = nearest(status.population, p; tol = 1e-16)
    η_ll = η_ll_max * rand()

    if d >= tol

        vv = (c_ll - v)
        # current-to-center
        y1 = y0 + (η_ll / norm(vv)) * vv
        y1 = correct(y1, bounds_ll)

        # approx
        r = Optim.optimize(
            z -> problem.f(p, z),
            problem.bounds_ll[1, :],
            problem.bounds_ll[2, :],
            y1,
            Optim.Fminbox(Optim.BFGS()),
        )
        f_calls += r.f_calls
        q = r.minimizer
        fpq = r.minimum
    else
        q = y0
        fpq = problem.f(p, q)
        f_calls += 1
    end

    return LLResult(q, fpq; f_calls = f_calls)

end


optimize(F, f, bounds_ul, bounds_ll, method::QBCA) =
    optimize(F, f, bounds_ul, bounds_ll, method)
