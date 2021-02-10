################################################################################
#                       L o w e r      L e v e l
################################################################################


########
#
# pensar el nivel de abajo como moead en problemas de un solo objetivo
#
#
######


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
        push!(B, I[2:T+1])
    end

    return B
end


function optimize_ll_population(X, problem, parameters, status, information, options)
    B = get_closest_vectors(X, parameters.T)

    a = problem.bounds_ll[1,:]
    b = problem.bounds_ll[2,:]
    D_ll = length(a)
    N = parameters.N

    Y = a' .+ (b - a)' .* rand(N, D_ll) 

    ll_population = [Metaheuristics.create_child(Y[i,:], problem.f(X[i],Y[i,:])) for i in 1:N]
    status.f_calls += N

    K = min(parameters.T, parameters.K)

    # main loop
    for gen in 1:options.ll_iterations
        #@show gen
        C = zeros(N, D_ll)
        ids_worst = Int[]

        # loop for each solution in ll population
        for i in 1:N

            if rand() < parameters.δ
                neighbors = shuffle(B[i])[1:K]
            else
                neighbors = shuffle(1:N)[1:K]
            end

            mass = Metaheuristics.fvals(ll_population[neighbors])
            # mass
            mass = 2maximum(abs.(mass)) .- mass
            mass = mass / sum(mass)

            # worst element in neighbors
            push!(ids_worst, neighbors[argmin(mass)])

            # center of mass
            C[i,:] = sum(Metaheuristics.positions(ll_population[neighbors]) .* mass, dims=1)
        end
        # display(C)

        Y_worst = Y[ids_worst, :]
        η = parameters.η_max*rand(N)

        Y_new = Y +  η .* (C - Y_worst)
        # fix solutions
        mask = .!(a' .<= Y_new .<= b')
        # println("mask ", sum(sum(mask, dims=2) .> 0 ) / N)
        
        Y_new[ mask ] = C[mask]

        #sleep(0.05)
        #display(Y)

        for i in 1:N
            new_sol = Metaheuristics.create_child(Y_new[i,:], problem.f(X[i], Y_new[i,:]))
            status.f_calls += 1
            j = ids_worst[i]
            if Metaheuristics.is_better_eca(new_sol, ll_population[j])
                ll_population[j] = new_sol
                Y[j,:] = Y_new[i,:]
                # println("gen: ", gen, "\ti: ", i, "\tf: ", new_sol.f)
            end
        end

        XX = hcat(X...)'
        plt.scatter(XX[:,end], Y[:,end], title="Gen: $gen", xlim=(-10, 10), ylim=(-10,10))
        #plt.scatter!(C[:,1], C[:,end], label="center")
        
        plt.gui()
    end

    return ll_population

    
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

    ll_pop = optimize_ll_population(X, problem, parameters, status, information, options)
    options.debug && display(Metaheuristics.positions(ll_pop))


    status.population = []
    for i in 1:length(X)
        x = X[i]
        y = ll_pop[i].x
        new_sol = generateChild(x, y, problem.F(x,y), problem.f(x, y))
        push!(status.population, new_sol)
    end

    status.stop = true


end


################################################################################
#                    I T E R A T I O N
################################################################################

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









mutable struct SABO2
    N::Int64
    K::Int64
    T::Int
    η_max::Float64
    λ::Float64
    δ::Float64
end

function SABO2(
    D::Int;
    K = 3,
    N = K * D,
    T = round(Int, 0.1N),
    η_max = 1.2,
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

    sabo = SABO2(N, K, T, η_max, λ, δ)

    algorithm = Algorithm(
        sabo;
        initialize! = initialize_SABO2!,
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
