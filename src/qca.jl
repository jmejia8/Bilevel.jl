#####################################################
#
#         BCA optimizer
#
#####################################################

function getU(P::Array, K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N-K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j .+ 1]
    end

    return P[U_ids]
end

function fitnessToMass(fitness::Array{Float64,2}, searchType::Symbol)
    m_ul = minimum(fitness[:,1])
    m_ll = minimum(fitness[:,2])
    
    m_ul < 0 && (fitness[:,1] = 2abs(m_ul) .+ fitness[:,1])
    m_ll < 0 && (fitness[:,1] = 2abs(m_ll) .+ fitness[:,2])
  

    if searchType == :minimize
        fitness[:,1] = 2maximum(fitness[:,1]) .- fitness[:,1]
        fitness[:,2] = 2maximum(fitness[:,2]) .- fitness[:,2]
    end

    return fitness
end


function getMass(U::Array, V::Array, α, β, searchType::Symbol)
    n = length(U)

    fitness = zeros(Float64, n, 2)
    
    for i = 1:n
        fitness[i, 1] = U[i].F + α*U[i].f
        fitness[i, 2] = V[i].f + α*U[i].F
    end

    return fitnessToMass(fitness, searchType)
end


function center(U::Array, V::Array, mass::Array{Float64,2})
    d_ul = length(U[1].x)
    d_ll = length(V[1].y)

    c_ul = zeros(Float64, d_ul)
    c_ll = zeros(Float64, d_ll)
    
    for i = 1:length(U)
        c_ul += mass[i,1] * U[i].x
        c_ll += mass[i,2] * V[i].y
    end

    return c_ul / sum(mass[:,1]), c_ll / sum(mass[:,2])
end

function center(U::Array, V::Array, α, β, searchType::Symbol)
    n = length(U)

    mass = getMass(U, V, α, β, searchType)

    a, b = center(U, V, mass)
    return a, b, getWorstInd(U), argmin(mass[:,2])
end

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

function optimize(F_ul::Function, # upper level objective function
                  f_ll::Function; # lower level objective function
                  bounds_ul::Array,
                  bounds_ll::Array,
                  stop_criteria::Function = x -> false,
                  search_type::Symbol = :minimize,

                  # general Options
                  options::Options = Options(),

                  # BCA parameters
                  k::Int = 3,
                  N::Int = 2k * size(bounds_ul, 2),
                  η_max::Real = 2.0,
                  max_evals_ul::Int = 1000size(bounds_ul, 2),
                  desired_accu = 1e-4,
                  
                  # # upper level restrictions
                  # G::Function  = (x, y) -> 0.0,
                  # H::Function  = (x, y) -> 0.0,

                  # # lower level restrictions
                  # g::Function  = (x, y) -> 0.0,
                  # h::Function  = (x, y) -> 0.0,

                  )

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

    convergence = [deepcopy(Population)]

    stop = false

    α = 0.05

    # start search
    while !stop
        I_ul = randperm(N)
        I_ll = randperm(N)

        success_rate = 0

        for i in 1:N

            # current
            x = Population[i].x
            y = Population[i].y

            # generate U masses
            U = getU(Population, k, I_ul, i, N)
            V = getU(Population, k, I_ll, i, N)
            
            # generate center of mass
            c_ul, c_ll, u_worst, v_worst = center(U, V, α, α, search_type)

            # stepsize
            η_ul = η_max * rand()
            η_ll = 0.5 * rand()

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
                r = Optim.optimize( z -> f(p, z), y1, Optim.BFGS())
                nevals_ll += r.f_calls
                q = r.minimizer
            else
                q = y0
            end
            
            q = correct(q, bounds_ll)

            sol = generateChild(p, q, F(p, q), f(p, q))
            nevals_ll += 1

            if sol ≺ Population[i]
                Population[getWorstInd(Population)] = sol
                success_rate += 1

                if sol ≺ best
                    best = sol
                    # push!(convergence, best)

                    stop = abs(best.f) < desired_accu && abs(best.F) < desired_accu
                    stop && break
                end
            end

        end


        push!(convergence, deepcopy(Population))

        iteration += 1

        stop = stop || (success_rate/N) < 0.01 || nevals_ul > max_evals_ul

    end


    return convergence, best, iteration, nevals_ul, nevals_ll
end
