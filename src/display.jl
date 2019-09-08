import Base.Multimedia.display

function display(solution::xf_indiv)
    @printf("| UL Fun: %g\n", solution.F)
    @printf("| LL fun: %g\n", solution.f)
    @printf("| UL Vec: "); display(solution.x)
    println("")
    @printf("| LL Vec: "); display(solution.y)
    println("")


end

function display(status::State)

    if typeof(status.best_sol) != xf_indiv
        return println(status)
    end

    println("+=========== STATE ==========+")
    @printf("| Iter.: %.0f\n", status.iteration)
    display(status.best_sol)
    

    # upper level parameters
    @printf("| UL FEs: %.0f\n", status.F_calls)

    # upper level parameters
    @printf("| LL FEs: %.0f\n", status.f_calls)
    @printf("|   Time: %0.4f secs.\n", status.final_time - status.initial_time)
    println("|   Stop: ", status.stop_msg)
    println("+============================+")
end

function display(method::QBCA)
    println("+=========== QBCA ===========+")
    @printf("|     k: %d\n", method.k)
    @printf("|     N: %d\n", method.N)
    @printf("|  η_ul: %g\n", method.η_ul_max)
    @printf("|  η_ll: %g\n", method.η_ll_max)
    @printf("| s_min: %g\n", method.s_min)

    # general Options
    @printf("|    iterations: %d\n", method.options.iterations)
    @printf("| F calls limit: %d\n", method.options.F_calls_limit)
    @printf("| f calls limit: %d\n", method.options.f_calls_limit)
    println("+============================+")
end

# show(status::State) = display(status)