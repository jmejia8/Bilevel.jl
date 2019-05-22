struct ParametersTemplate
    N::Int # population size
end


function initialize!(problem::Problem, parameters, status::State)
    # some stuff here
end

function update_state!(problem::Problem,status::State,information::Information,options::Options,t::Int)
    # some stuff here
end

function lower_level_optimizer(problem::Problem,status::State,information::Information,options::Options,t::Int)
    # some stuff here
end

function is_better(solution_1, solution_2) # solution_1 is better that solution_2 
    solution_1 ≺ solution_2
end

function stop_criteria(status::State, information::Information, options::Options)
    stop_check(status, information, options)
end

function final_stage!(status::State, information::Information, options::Options)
    # some stuff here
end
