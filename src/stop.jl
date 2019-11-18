function ull_call_limit_stop_check(status::State, information::Information, options::Options)
    status.stop =  status.F_calls >= options.F_calls_limit || status.f_calls >= options.f_calls_limit
    
    if status.stop
        status.stop_msg = "ul_ll_call_limit"
        options.debug && @info("Stopped since ull_call_limit was met.")
    end
    
    status.stop
end

function iteration_stop_check(status::State, information::Information, options::Options)
    status.stop =  status.iteration >= options.iterations
    
    if status.stop
        status.stop_msg = "iteration"
        options.debug && @info("Stopped since iteration limit was met.")
    end
    
    status.stop
end

function ul_accuracy_stop_check(status::State, information::Information, options::Options)
    status.stop =  information.F_optimum != NaN && abs(status.best_sol.F - information.F_optimum) < options.F_tol
    
    if status.stop
        status.stop_msg = "ul_accuracy"
        options.debug && @info("Stopped since ul_accuracy was met.")
    end
    
    status.stop
end

function ll_accuracy_stop_check(status::State, information::Information, options::Options)
    status.stop = information.f_optimum != NaN && abs(status.best_sol.f - information.f_optimum) < options.f_tol
    
    if status.stop
        status.stop_msg = "ll_accuracy"
        options.debug && @info("Stopped since ll_accuracy was met.")
    end
    
    status.stop
end

function ul_varF_stop_check(status::State, information::Information, options::Options)
    status.stop =  var(map( s->s.F, status.population )) ≈ 0.0
    
    if status.stop
        status.stop_msg = "ul_varF"
        options.debug && @info("Stopped since ul_varF was met.")
    end
    
    status.stop
end

function ul_diversity_stop_check(status::State, information::Information, options::Options)
    status.stop =  var(map( s->s.F, status.population )) ≈ 0.0
    
    if status.stop
        status.stop_msg = "ul_diversity"
        options.debug && @info("Stopped since ul_varF was met.")
    end
    
    status.stop
end


function accuracy_stop_check(status::State, information::Information, options::Options)
    status.stop = ul_accuracy_stop_check(status, information, options)  && ll_accuracy_stop_check(status, information, options)

    if status.stop
        status.stop_msg = "ul_ll_accuracy"
    else
        status.stop_msg = ""
    end

    status.stop
end


function stop_check(status::State, information::Information, options::Options)
    ull_call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)    ||
    ul_varF_stop_check(status, information, options)    ||
    accuracy_stop_check(status, information, options)
end