module BCAOperators

import  ..stop_check, ..Selection


function initialize!(problem, parameters)
    # some stuff here
end

function update_state!(problem,status,information,options,t)
    # some stuff here
end

function lower_level_optimizer(x,problem,status,information,options,t)
    # some stuff here
end

function is_better(solution_1, solution_2) # solution_1 is better that solution_2 
    solution_1 ≺ solution_2
end

function stop_criteria(status, information, options)
    stop_check(status, information, options)
end

function final_stage!(status, information, options)
    # some stuff here
end

end