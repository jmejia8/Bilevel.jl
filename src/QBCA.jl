function nearest(P, x; tol = 1e-5)
    x_nearest = P[1].x
    y = P[1].y
    d = Inf

    for sol = P
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

function optimize(parameters::QBCA)
    F_ul      = parameters.F
    f_ll      = parameters.f
    bounds_ul = parameters.bounds_ul
    bounds_ll = parameters.bounds_ll
    stop_criteria = parameters.stop_criteria
    search_type   = parameters.search_type

    options = parameters.options

    k = parameters.k
    N = parameters.N
    η_ul_max = parameters.η_ul_max
    η_ll_max = parameters.η_ll_max

    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 
    nevals_ul = 0

    F(x, y) = begin
        nevals_ul += 1
        F_ul(x, y)
    end

    nevals_ll = 0
    nevals_lll = 0

    f(x, y) = begin
        nevals_lll += 1
        f_ll(x, y)
    end

    # general parameters
    D = D_ul + D_ll
    # N = D < 5 ? 5*D : k*D

    # tolerance for y proximity
    tol = 1e-16

    # initialize population
    Population, f_calls = init_population(F, f, N, bounds_ul, bounds_ll)

    nevals_ll += f_calls

    # current generation
    iteration = 0

    # best solution
    best = getBest(Population)
    status = State(best)

    status.F_calls   = nevals_ul
    status.f_calls   = nevals_lll
    status.iteration = iteration

    convergence = [deepcopy(Population)]

    stop = false

    α = 0.05
    β = 0.05

    # start search
    while !stop
        I_ul = randperm(N)
        I_ll = randperm(N)

        status.success_rate = 0

        for i in 1:N

            # current
            x = Population[i].x
            y = Population[i].y

            # generate U masses
            U = getU(Population, k, I_ul, i, N)
            V = getU(Population, k, I_ll, i, N)
            
            # generate center of mass
            c_ul, c_ll, u_worst, v_worst = center(U, V, α, β, search_type)

            # stepsize
            η_ul = η_ul_max * rand()
            η_ll = η_ll_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            v = V[v_worst].y
            
            # current-to-center
            p = x + η_ul * (c_ul - u)
            p = correct(p, bounds_ul)

            y0, d = nearest(Population, p; tol = tol)
            
            if d >= tol

                vv = (c_ll - v)
                # current-to-center
                y1 = y0 + (η_ll/ norm(vv)) * vv
                y1 = correct(y1, bounds_ll)
                
                # approx
                r = Optim.optimize( z -> f(p, z), bounds_ll[1,:], bounds_ll[2,:], y1, Optim.Fminbox(Optim.BFGS()))
                status.f_calls += r.f_calls
                q = r.minimizer
                fpq = r.minimum
            else
                q = y0
                fpq = f(p, q)
                status.f_calls += 1
            end
            
            # q = correct(q, bounds_ll)

            sol = generateChild(p, q, F(p, q), fpq)
            status.F_calls += 1

            if sol ≺ Population[i]
                Population[getWorstInd(Population)] = sol
                status.success_rate += 1

                if sol ≺ best
                    status.best_sol = sol
                    # push!(convergence, status.best_sol)

                    stop = abs(status.best_sol.f) < options.f_tol && abs(status.best_sol.F) < options.F_tol
                    stop && break
                end
            end

        end


        push!(convergence, deepcopy(Population))

        status.iteration += 1

        stop = stop || (status.success_rate/N) < parameters.s_min || nevals_ul > parameters.F_calls_limit

    end


    return status#convergence, best, iteration, nevals_ul, nevals_ll
end
