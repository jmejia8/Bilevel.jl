module BCAOperators

import  ..stop_check, ..Selection, ..generateChild, ..LLResult, ..getU, ..Metaheuristics
import Random: randperm
include("utils.jl")

struct BCAFW
    N::Int64
    K::Int64
    η_max::Float64
end

BCAFW(;N=100, K = 7, η_max = 2.0) = BCAFW(N, K, η_max)

function initialize!(problem,parameters,status,information,options)
    p = problem

    a = p.bounds_ul[1,:]
    b = p.bounds_ul[2,:]

    X = a'  .+ (b - a)' .* rand(parameters.N, size(p.bounds_ul, 2))

    # population array
    population = []

    for i in 1:parameters.N
        x = X[i,:]
        ll_result = lower_level_optimizer(x, p, status, information, options, 0)

        y = ll_result.y

        child = generateChild(x, y, p.F(x, y), p.f(x, y))
        push!(population, child)

        status.f_calls += ll_result.f_calls
        status.F_calls += 1
    end

    individual = typeof(population[1])
    population = Array{individual, 1}(population)

    status.population = population
    status.best_sol = population[1]
end

function update_state!(problem,engine,parameters,status,information,options,t)
    I = randperm(parameters.N)
    c = zeros(size(problem.bounds_ul,2))


    is_better = engine.is_better

    i::Int = 1
    for sol in status.population
        x = sol.x

        ########################################################################
        # center of mass
        ########################################################################
        # generate U masses
        U = getU(status.population, parameters.K, I, i, parameters.N); i+=1
        
        # generate center of mass
        center!(c,U)

        ########################################################################
        # New upper level vector
        ########################################################################
        # stepsize
        η_ul = parameters.η_max * rand()

        # u: worst element in U
        u = worst(U, is_better)

        p = x + η_ul * (c - u)
        correctSol!(p, problem.bounds_ul[1,:], problem.bounds_ul[2,:])

        ########################################################################
        # Nested optimization solution
        ########################################################################
        ll_result = engine.lower_level_optimizer(p,problem,status,information,options,0)
        status.f_calls += ll_result.f_calls
        q = ll_result.y
        ########################################################################

        new_sol = generateChild(p, q, problem.F(p, q), ll_result.f)
        status.F_calls += 1
        ########################################################################
        # comparing new solution
        ########################################################################

        if is_better(new_sol, sol)
            status.population[getWorstInd(status.population, is_better)] = new_sol

            if is_better(new_sol, status.best_sol)
                status.best_sol = new_sol
            end
        end

        
        ########################################################################
        fill!(c, 0.0)
    end
end

function lower_level_optimizer(x,problem,status,information,options,t)
    D = size(problem.bounds_ll, 2)
    y, fy = Metaheuristics.eca( z -> problem.f(x, z), D;
                                limits=problem.bounds_ll,
                                K = 3, N = 3D,
                                showResults=false,
                                p_bin = 0,
                                p_exploit = 2,
                                canResizePop=false, 
                                max_evals=1000D)
    return LLResult(y, fy; f_calls = 1000D)
end

function is_better(solution_1, solution_2) # solution_1 is better that solution_2 
    Selection(solution_2, solution_1)
end

function stop_criteria(status, information, options)
    stop_check(status, information, options)
end

function final_stage!(status, information, options)
    # some stuff here
end

end