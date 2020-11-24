g_te(fxy, λ, z) = maximum( λ .* abs.(fxy - z) )

function update_z_ideal!(z, fx::Vector{Float64})
    for j in 1:length(z)
        if z[j] > fx[j]
            z[j] = fx[j]
        end
    end
end

function update_z_ideal!(z, ll_population)
    for i in 1:length(ll_population)
        update_z_ideal!(z, ll_population[i].f)
    end

end

function update_neighbors_ll!(X, parameters)
    N = length(X)
    d = zeros(N, N)

    empty!(parameters.B_ll)

    for i in 1:N
        for j in (i+1):N
            d[i, j] = norm(X[i] - X[j])
            d[j, i] = d[i, j]
        end

        I = sortperm(d[i,:])
        push!(parameters.B_ll, I[2:(parameters.T_ll+1)])
    end

    parameters.B_ll
end

function getU(B_ll::Array{Int}, K)
    I = shuffle(B_ll)
    return view(I, 1:K)
end

function get_mass(U, G_te)
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        v = U[i].sum_violations
        if v > 0.0
            fitness[i] = v
        else
            fitness[i] = G_te[i]
        end
    end

    return Metaheuristics.fitnessToMass(fitness, :minimize)

end

function lower_level_iteration!(X, ll_population, G_te, problem, status, information, options, t, parameters)
    D_ll = size(problem.bounds_ll,2)
    D_ul = size(problem.bounds_ul,2)

    a = problem.bounds_ll[1,:]
    b = problem.bounds_ll[2,:]

    # generate center of mass
    for i = 1:parameters.N
        x = X[i]
        y = ll_population[i].x
        sol = ll_population[i]
        # generate U masses (neighbors of y_i)
        U_ids = getU(parameters.B_ll[i], parameters.K_ll)
        U = view(ll_population, U_ids)

        mass = get_mass(U, G_te)
        c = Metaheuristics.center(U, mass)
        u_worst = argmin(mass)


        # stepsize
        η = parameters.η_max_ll * rand()

        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        y_new = y .+ η .* (c .- u)

        y_new = Metaheuristics.correct(y_new, c, a, b)

        # fxy, gxy, hxy = problem.f(x, y_new)

        sol_new = Metaheuristics.generateChild(y_new, problem.f(x, y_new))
        sol_new.sum_violations = Metaheuristics.violationsSum(sol_new.g, sol_new.h)
        status.f_calls += 1
        update_z_ideal!(parameters.w, sol_new.f)

        vio_new = sol_new.sum_violations
        vio_old = sol.sum_violations

        for j in parameters.B_ll[i]
            gg = g_te(sol_new.f, parameters.λ_ll[j], parameters.w)
            gg2 = g_te(sol.f, parameters.λ_ll[j], parameters.w)



            is_better_than_current = (vio_new < vio_old) || (vio_new == vio_old && gg <= gg2)

            if !is_better_than_current
                continue
            end

            ll_population[j] = sol_new
            G_te[j] = gg


        end


        # status.stop = engine.stop_criteria(status, information, options)
        status.stop && break
    end

end



function lower_level_optimizer_MBO(X, problem, status, information, options, t, parameters=nothing)
    if isnothing(parameters)
        error("Lower level optimizer require parameters")
    end

    # initialization
    a = problem.bounds_ll[1,:]
    b = problem.bounds_ll[2,:]


    ll_population = Metaheuristics.initializePop(y -> ([0.0, 0], [0.0], [0.0]), parameters.N, length(a), a, b)
    ff(f_ll, λ, z) = begin
        fy = f_ll[1]
        gy = f_ll[2]
        hy = f_ll[3]
        return (g_te(fy, λ, z), gy, hy)
    end

    for i in 1:length(ll_population)
        if isempty(parameters.λ_ll[i])
            parameters.λ_ll[i] = rand(range(0, 1, length=parameters.N), 2)
            parameters.λ_ll[i] /= sum(parameters.λ_ll[i])
        end

        if isempty(parameters.w)
            parameters.w = ones(2)*Inf
        end
        @show parameters.λ_ll[i]
        @show parameters.w

        eca = Metaheuristics.ECA()
        eca.options.f_calls_limit = 1000
        r = Metaheuristics.optimize(y-> ff(problem.f(X[i], y), parameters.λ_ll[i], parameters.w), problem.bounds_ll, eca)

        # ll_population[i].x = parameters.Ψ(X[i])
        ll_population[i].x = r.best_sol.x
        fxy, gxy, hxy = problem.f(X[i], ll_population[i].x)
        status.f_calls += r.f_calls + 1
        ll_population[i].f = fxy
        ll_population[i].g = gxy
        ll_population[i].h = hxy



        update_z_ideal!(parameters.w, fxy)
    end






    G_te = [ g_te(ll_population[i].f, parameters.λ_ll[i], parameters.w) for i in 1:parameters.N ]

    # find neighbors of each X[i]
    update_neighbors_ll!(X, parameters)

    # main loop
    iterations = 0
    for t in 1:iterations
        lower_level_iteration!(X, ll_population, G_te, problem, status, information, options, t, parameters)
    end

    # x = map(s -> s.f[1], ll_population)
    # y = map(s -> s.f[2], ll_population)
    # plt = scatterplot(x, y)
    # display(plt)

    return Bilevel.LLResult(ll_population, nothing; f_calls = 0)


end
