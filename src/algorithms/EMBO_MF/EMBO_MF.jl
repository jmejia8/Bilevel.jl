include("lower-level.jl")
mutable struct MBO
    upper_level_parameters
    lower_level_parameters
    n_followers::Int
    Ψ::Function
end


function MBO(;
        upper_level_parameters = nothing,
        lower_level_parameters = nothing,
        Ψ = x -> x,
        n_followers = 1,
        options = Options(),
        information = Information()
)



    parameters = MBO(upper_level_parameters, lower_level_parameters, n_followers, Ψ)

    algorithm = Bilevel.Algorithm(
        parameters;
        initialize! = initialize_MBO!,
        update_state! = update_state_MBO!,
        lower_level_optimizer = lower_level_optimizer_MBO,
        is_better = is_better_MBO,
        stop_criteria = stop_criteria_MBO,
        final_stage! = final_stage_MBO!,
        options = options,
        information = information,
    )



    algorithm
end

function ll_decision_based_on_nadir(F, x, followers_result_)
    followers_result = followers_result_.y
    n_followers = length(followers_result)
    Y = []
    for i in 1:n_followers
        @show i
        display(followers_result[i])
        pareto_solutions = (followers_result[i].best_sol) 

        # get the ideal
        z = fill(Inf, length(pareto_solutions[1].f))
        Metaheuristics.update_reference_point!(z, pareto_solutions)

        j = argmin( [ norm( z - sol.f ) for sol in pareto_solutions ] )

        push!(Y, Metaheuristics.get_position(pareto_solutions[j]))
    end

    return Y
end


function ll_decision_maker(F, x, Y)
    n_followers = length(Y)
    Fx = []
    # for i in 1:n_followers
    #     no_dominated = Y[i] # for the i-th follower
    #     for j in 1:length(no_dominated)
    #         F(x,)
    #     end
    # end
end

function initialize_MBO!(problem, engine, parameters, status, information, options)

    D_ul = size(problem.bounds_ul, 2)
    N = parameters.upper_level_parameters.parameters.N
    a = problem.bounds_ul[1,:]
    b = problem.bounds_ul[2,:]
    X = [ a + (b - a).*rand(D_ul) for i in 1:N ] 

    status.population = []

    for x in X
        options.debug && @show x
        ll_results = lower_level_optimizer_MBO(x, problem, status, information, options, 0;
                                        parameters.lower_level_parameters,
                                        parameters.n_followers )

        Y = ll_decision_based_on_nadir(problem.F, x, ll_results)
        child = generateChild(x, Y, problem.F(x, Y), (zeros(0), [0.0], [0.0] ))
        child.other = ll_results

        push!(status.population, child)
        break
    end

    status.stop = true
  
    status.final_time = time()
end



function update_state_MBO!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        t,
       )



end


function is_better_MBO(
    New::xFgh_indiv,
    Old::xFgh_indiv;
    ε = 0.0,
)

    old_vio = Metaheuristics.violationsSum(Old.G, Old.H, ε = ε) # upper level
    old_vio += Metaheuristics.violationsSum(Old.g, Old.h, ε = ε) # lower level

    new_vio = Metaheuristics.violationsSum(New.G, New.H, ε = ε) # upper level
    new_vio += Metaheuristics.violationsSum(New.g, New.h, ε = ε) # lower level

    if new_vio < old_vio
        return true
    elseif new_vio > old_vio
        return false
    end

    for i in 1:length(Old.F)
        if Old.F[i] < New.F[i]
            return false
        end
    end

    return true
end


function stop_criteria_MBO(status::State, information::Information, options::Options)
    Bilevel.ull_call_limit_stop_check(status, information, options) ||
    Bilevel.iteration_stop_check(status, information, options)
end

function final_stage_MBO!(status, information, options)
    status.final_time = time()

end


