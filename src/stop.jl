function ull_call_limit_stop_check(status::State, information::Information, options::Options)
    cond =  status.F_calls >= options.F_calls_limit || status.f_calls >= options.f_calls_limit
    options.debug && cond && @info("Stopped since ull_call_limit was met.")
    cond
end

function iteration_stop_check(status::State, information::Information, options::Options)
    cond =  status.iteration >= options.iterations
    options.debug && cond && @info("Stopped since iteration limit was met.")
    cond
end

function ul_accuracy_stop_check(status::State, information::Information, options::Options)
    cond =  information.F_optimum != NaN && abs(status.best_sol.F - information.F_optimum) < options.F_tol
    options.debug && cond && @info("Stopped since ul_accuracy was met.")
    cond
end

function ll_accuracy_stop_check(status::State, information::Information, options::Options)
    information.f_optimum != NaN && abs(status.best_sol.f - information.f_optimum) < options.f_tol
end

function ul_varF_stop_check(status::State, information::Information, options::Options)
    cond =  var(map( s->s.F, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since ul_varF was met.")
    cond
end

function ul_diversity_stop_check(status::State, information::Information, options::Options)
    cond =  var(map( s->s.F, status.population )) ≈ 0.0
    options.debug && cond && @info("Stopped since ul_varF was met.")
    cond
end



function stop_check(status::State, information::Information, options::Options)
    ull_call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)    ||
    ul_varF_stop_check(status, information, options)    ||
    (ul_accuracy_stop_check(status, information, options)  && ll_accuracy_stop_check(status, information, options))
end