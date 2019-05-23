module BCAOperators

import  ..stop_check, ..Selection

function initialize!(problem, parameters, status)
    p = problem
    population, f_calls = init_population(p.F, p.f, parameters.N, p.bounds_ul,p.bounds_ll)
    
    status.population = population
    status.f_calls = f_calls

end

function update_state!(problem,status,information,options,t)
    # some stuff here
end

function lower_level_optimizer(problem,status,information,options,t)
    # some stuff here
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