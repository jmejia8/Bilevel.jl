function ull_call_limit_stop_check(status::State, information::Information, options::Options)
    status.F_calls >= options.F_calls_limit || status.f_calls >= options.f_calls_limit
end

function iteration_stop_check(status::State, information::Information, options::Options)
    status.iteration >= options.iterations
end

function ul_accuracy_stop_check(status::State, information::Information, options::Options)
    information.F_optimum != NaN && abs(status.best_sol.F - information.F_optimum) < options.F_tol
end

function ll_accuracy_stop_check(status::State, information::Information, options::Options)
    information.f_optimum != NaN && abs(status.best_sol.f - information.f_optimum) < options.f_tol
end

function stop_check(status::State, information::Information, options::Options)
    ull_call_limit_stop_check(status, information, options) ||
    iteration_stop_check(status, information, options)    ||
    (ul_accuracy_stop_check(status, information, options)  && ll_accuracy_stop_check(status, information, options))
end