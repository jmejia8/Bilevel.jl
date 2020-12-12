################################################################################
#                       L o w e r      L e v e l
################################################################################




################################################################################
#                    I n i t i l i a z a t i o n
################################################################################

import Plots
plt = Plots

function get_closest_vectors(X, T)
    N = length(X)
    distances = zeros(N, N)
    λ = X
    B = []
    for i in 1:N
        for j in (i+1):N
            distances[i, j] = norm(λ[i] - λ[j])
            distances[j, i] = distances[i, j]
        end
        I = sortperm(distances[i, :])
        if T > 1
            push!(B, I[2:T+1])
        else 
            push!(B, I[2])
        end
    end

    return B
end


function optimize_ll_population(X, problem, parameters, status, information, options; Y = [])
    B = []

    a = problem.bounds_ll[1,:]
    b = problem.bounds_ll[2,:]
    D = length(a)
    N = parameters.N

    if isempty(Y) || isnothing(Y)
        Y = [ a + (b - a) .* rand(D) for i in 1:N]
    end


    population = [Metaheuristics.create_child(Y[i], problem.f(X[i],Y[i])) for i in 1:N]
    status.f_calls += N

    K = min(parameters.T, parameters.K)
    
    # main loop
    tt = 0
    for gen in 1:options.ll_iterations
        tt += 1
        n_improvemts = 0
        # loop for each solution in ll population
        for i in 1:N

            P_idx = collect(1:N)

            # select participats
            r1 = i
            r2 = rand(P_idx)
            while r1 == r2
                r2 = rand(P_idx)
            end

            r3 = rand(P_idx)
            while r3 == r1 || r3 == r2
                r3 = rand(P_idx)
            end

            a = population[r1].x
            b = population[r2].x
            c = population[r3].x


            # binomial crossover
            v = zeros(D)
            j_rand = rand(1:D)

            # binomial crossover
            for j = 1:D
                # binomial crossover
                if rand() < parameters.CR
                    v[j] = a[j] + parameters.η_max * (b[j] - c[j])
                else
                    v[j] = a[j]
                end
                # polynomial mutation

                if rand() < parameters.p_m
                    r = rand()
                    if r < 0.5
                        σ_k = (2.0 * r)^(1.0 / (parameters.η + 1)) - 1
                    else
                        σ_k = 1 - (2.0 - 2.0 * r)^(1.0 / (parameters.η + 1))
                    end
                    v[j] = v[j] + σ_k * (problem.bounds_ll[2,j] - problem.bounds_ll[1,j])
                end
            end

            v = Metaheuristics.replace_with_random_in_bounds!(v, problem.bounds_ll)

            # instance child
            h = Metaheuristics.generateChild(v, problem.f(X[i], v))
            status.f_calls += 1

            if Metaheuristics.is_better_eca(h, population[i])
                population[i] = h
                n_improvemts += 1
            end



        end

        if n_improvemts == 0
            options.debug && @info "Early stop at ll $tt"
            break
        end


    end

    return population

    
end




function initialize_SABO2!(
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
    D_ul = length(a)

    X = [ a + (b - a) .* rand(D_ul) for i in 1:parameters.N]

    status.population = []

    for i in 1:length(X)
        x = X[i]
        res = Metaheuristics.optimize(yy -> problem.f(x,yy), problem.bounds_ll, Metaheuristics.ECA( ;options=Metaheuristics.Options(f_calls_limit=7000) ) )
        println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        println(">>>>>>>>>>>>>>>>>>>>>>>>  ",i,"  >>>>>>>>>>>>>>>>>>>>>>>")
        display(res)
        println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        y = Metaheuristics.minimizer( res )
        new_sol = generateChild(x, y, problem.F(x,y), problem.f(x, y))
        status.f_calls += res.f_calls
        status.F_calls += 1
        push!(status.population, new_sol)
        if isnothing(status.best_sol) || is_better_SABO(new_sol, status.best_sol)
            status.best_sol = new_sol
        end
    end




end


################################################################################
#                    I T E R A T I O N
################################################################################

function update_state_SABO2_simple!(
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
    c = zeros(size(problem.bounds_ul,2))


    i::Int = 1
    X = []
    for sol in status.population
        x = sol.x

        ########################################################################
        # center of mass
        ########################################################################
        # generate U masses
        U = getU(status.population, parameters.K, I, i, N); i+=1
        
        # generate center of mass
        c = BCAOperators.center!(c,U)

        ########################################################################
        # New upper level vector
        ########################################################################
        # stepsize
        η_ul = parameters.η_max * rand()

        # u: worst element in U
        u = BCAOperators.worst(U, engine.is_better)

        p = x + η_ul * (c - u)
        BCAOperators.correctSol!(p, c, problem.bounds_ul[1,:], problem.bounds_ul[2,:])

        push!(X, p)

    end

    B = get_closest_vectors(X,1)
    Y = map(s -> s.y, status.population[B])
    mask = rand(length(X)) .< 0.3
    Y[mask] = [ problem.bounds_ll[1,:] + (problem.bounds_ll[2,:] - problem.bounds_ll[1,:] ) .* rand(length(Y[1])) for i in 1:sum(mask) ]
    ll_pop = optimize_ll_population(X, problem, parameters, status, information, options;Y)

    for i in 1:N
        x = X[i]
        y = ll_pop[i].x

        new_sol = generateChild(x, y, problem.F(x,y), problem.f(x, y))
        status.f_calls +=1
        status.F_calls += 1

        if is_better_SABO(new_sol, status.population[i])
            status.population[i] = new_sol
        end
    end
    status.stop = engine.stop_criteria(status, information, options)
    
end



function update_state_SABO2!(
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

    update_state_SABO2_simple!(
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









mutable struct SABO2
    N::Int64
    K::Int64
    T::Int
    η_max::Float64
    CR::Float64
    p_m::Float64
    η::Float64
    λ::Float64
    δ::Float64
end

function SABO2(
    D::Int;
    K = 3,
    N = K * D,
    T = 20,
    η_max = 0.5,
    CR = 1.0,
    p_m = 1/5,
    η = 20.0,
    λ = 1e-5,
    δ = 0.9,
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
            ll_iterations = 500,
            F_tol = 1e-2,
            f_tol = 1e-3,
            store_convergence = true,
        )
    end

    # if isnothing(information)
    #     information = Information(F_optimum = 0.0, f_optimum = 0.0)
    # end

    sabo = SABO2(N, K, T, η_max, CR, p_m, η, λ, δ)

    algorithm = Algorithm(
        sabo;
        initialize! = initialize_SABO!,
        update_state! = update_state_SABO2!,
        lower_level_optimizer = lower_level_optimizer_SABO,
        is_better = is_better_SABO,
        stop_criteria = stop_criteria_SABO,
        final_stage! = final_stage_SABO!,
        options = options,
        information = information,
    )



    algorithm

end

optimize(F, f, bounds_ul, bounds_ll, method::SABO2) =
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
#     method = SABO2Module.SABO2(D_ul; K = K,
#                                      N = N,
#                                      η_max = η_max,
#                                      options = options,
#                                      information = information)
# end
