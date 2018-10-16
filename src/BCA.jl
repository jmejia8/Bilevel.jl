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
        fitness[i, 2] = β*V[i].F + V[i].f
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
    return a, b, getWorstInd(U, searchType), argmin(mass[:,2])
end

################################################################################
################################################################################
################################################################################

function optimize(F::Function, # upper level objective function
                  f::Function; # lower level objective function
                  bounds_ul::Matrix,
                  bounds_ll::Matrix,
                  stop_criteria::Function = x -> false,
                  search_type::Symbol = :minimize,

                  # general Options
                  options::Options = Options(),

                  # BCA parameters
                  k::Int = 3,
                  N::Int = k * size(bounds_ul, 2) * size(bounds_ll, 2),
                  η_max::Int = 2.0,
                  
                  # # upper level restrictions
                  # G::Function  = (x, y) -> 0.0,
                  # H::Function  = (x, y) -> 0.0,

                  # # lower level restrictions
                  # g::Function  = (x, y) -> 0.0,
                  # h::Function  = (x, y) -> 0.0,

                  )

    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 

    # general parameters
    D = D_ul + D_ll
    N = D < 5 ? N : κ*D

    # initialize population
    Population = init_population(F, f, N, bounds_ul, bounds_ll)

    # current generation
    iteration = 0

    # best solution
    best = getBest(Population, search_type)

    convergence = [best]

    # start search
    for iteration = 1:100
        I_ul = randperm(N)
        I_ll = randperm(N)

        for i in 1:N

            # current
            x = Population[i].x
            y = Population[i].y

            # generate U masses
            U = getU(Population, κ, I_ul, i, N)
            V = getU(Population, κ, I_ll, i, N)
            
            # generate center of mass
            c_ul, c_ll, u_worst, v_worst = center(U, V, α, β, search_type)

            # stepsize
            η_ul = η_max * rand()
            η_ll = η_max * rand()

            # u: worst element in U
            u = U[u_worst].x
            v = V[v_worst].y
            
            # current-to-center
            p = x + η_ul * (c_ul - u)
            p = correct(p, bounds_ul)
            
            r = optimize( z -> f(p, z), y, BFGS())
            q = r.minimizer

            sol = generateChild(p, q, F(p, q), f(p, q))

            # replace worst element
            if sol ≺ Population[i]
                Population[getWorstInd(Population, search_type)] = sol

                if is_better_mass(sol, best, search_type)
                    best = sol
                    push!(convergence, best)
                end
            end

        end

    end



    if true
        println("+----------------------------------+")
        println("|          HBO results             |")
        println("+----------------------------------+")
        printResults(best, Population, t, nevals_ul, nevals_ll)
        println("+----------------------------------+")
    end

    return best.x, best.y, best, nevals_ul, Population, convergence
end