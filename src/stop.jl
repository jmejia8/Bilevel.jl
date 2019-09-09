function ull_call_limit_stop_check(status::State, information::Information, options::Options)
    cond =  status.F_calls >= options.F_calls_limit || status.f_calls >= options.f_calls_limit
    options.debug && cond && @info("Stopped since ull_call_limit was met.")
    status.stop = cond
    status.stop_msg = "ul_ll_call_limit"
    cond
end

function iteration_stop_check(status::State, information::Information, options::Options)
    cond =  status.iteration >= options.iterations
    options.debug && cond && @info("Stopped since iteration limit was met.")
    status.stop = cond
    status.stop_msg = "iteration"
    cond
end

function ul_accuracy_stop_check(status::State, information::Information, options::Options)
    cond =  information.F_optimum != NaN && abs(status.best_sol.F - information.F_optimum) < options.F_tol
    options.debug && cond && @info("Stopped since ul_accuracy was met.")
    status.stop = cond
    status.stop_msg = "ul_accuracy"
    cond
end

function ll_accuracy_stop_check(status::State, information::Information, options::Options)
    cond = information.f_optimum != NaN && abs(status.best_sol.f - information.f_optimum) < options.f_tol
    status.stop = cond
    status.stop_msg = "ll_accuracy"
    cond
end

function ul_varF_stop_check(status::State, information::Information, options::Options)
    cond =  var(map( s->s.F, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since ul_varF was met.")
    status.stop = cond
    status.stop_msg = "ul_varF"
    cond
end

function ul_diversity_stop_check(status::State, information::Information, options::Options)
    cond =  var(map( s->s.F, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since ul_varF was met.")
    status.stop = cond
    status.stop_msg = "ul_diversity"
    cond
end



function stop_check(status::State, information::Information, options::Options)
    ull_call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)    ||
    ul_varF_stop_check(status, information, options)    ||
    (ul_accuracy_stop_check(status, information, options)  && ll_accuracy_stop_check(status, information, options))
end