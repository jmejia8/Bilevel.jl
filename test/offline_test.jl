using BilevelBenchmark

include("../src/externals.jl")
include("../src/structures.jl")
include("../src/display.jl")
include("../src/tools.jl")
include("../src/operators.jl")


include("../src/BCA.jl")
include("../src/QBCA.jl")

using CSVanalyzer
import DelimitedFiles.writedlm


function getBilevel(fnum)

    if fnum == 0
        L(w, a, b) = a + w^2 - a*cos( π * b * w)
        f0(x,y; n=length(y)-1) = abs( 4x[n+1] - x[n+2] ) * y[end]^2 + sum(L.( abs.(x[1:n]-y[1:n]), x[n+1], x[n+2]))
        F0(x,y; n=length(y)-1) = abs( 4x[n+1] - x[n+2] ) / (1+y[end]^2) + sum(L.( abs.(x[1:n]-y[1:n]), x[n+1], x[n+2]))

        n = 3
        bounds_ul = zeros(2, n+2)
        bounds_ul[1,1:n] = -10ones(n)
        bounds_ul[2,1:n] = 10ones(n)
        bounds_ul[1, n+1:n+2] = [2, 2.0]
        bounds_ul[2, n+1:n+2] = [10, 10.0]
        
        bounds_ll = ones(2, n+1)
        bounds_ll[1,:] = -10ones(n+1)
        bounds_ll[2,:] = 10ones(n+1)

        return F0, f0, bounds_ul, bounds_ll

    end

    D_ul = D_ll = 5

    bounds_ul, bounds_ll = bilevel_ranges(D_ul, D_ll, fnum)

    # leader
    F(x::Array{Float64}, y::Array{Float64}) =  bilevel_leader(x, y, fnum)

    # follower
    f(x::Array{Float64}, y::Array{Float64}) = bilevel_follower(x, y, fnum)

    return F, f, bounds_ul, bounds_ll

end

function test()



    for fnum = 1:8
        F, f, bounds_ul, bounds_ll = getBilevel(fnum)
        method = QBCA(size(bounds_ul, 2); options = Options(F_tol=1e-4, f_tol=1e-4))
        result = optimize(F, f, bounds_ul, bounds_ll, method)
        b   = result.best_sol
        iters  = result.iteration
        nevals_ul = result.F_calls
        nevals_ll = result.f_calls
        

        @printf("SMD%d \t F = %e \t f = %e \t ev_ul = %d \t ev_ll = %d\n", fnum, b.F, b.f, nevals_ul, nevals_ll)
    end



end

function configure(output_dir="./")
    all_data_dir = string(output_dir, "output")
    summary_data_dir = string(output_dir, "summary")
    
    !isdir(all_data_dir)     && mkdir(all_data_dir)
    !isdir(summary_data_dir) && mkdir(summary_data_dir)

    # println("Using ", Threads.nthreads(), " threads.")
end

function saveData(output_dir="./")

    configure()

    NFUNS = 8
    NRUNS = 31

    errors_UL = zeros(NFUNS, NRUNS)
    errors_LL = zeros(NFUNS, NRUNS)

    evals_UL = zeros(NFUNS, NRUNS)
    evals_LL = zeros(NFUNS, NRUNS)


    for fnum = 1:8
        
        F, f, bounds_ul, bounds_ll = getBilevel(fnum)
        for i = 1:NRUNS

            P, best, iters, nevals_ul, nevals_ll = optimize(F, f, bounds_ul = bounds_ul, bounds_ll=bounds_ll)
        
            @printf("SMD%d \t F = %e \t f = %e \t ev_ul = %d \t ev_ll = %d\n", fnum, best.F, best.f, nevals_ul, nevals_ll)

            errors_UL[fnum, i] = best.F; errors_LL[fnum, i] = best.f
            evals_UL[fnum, i] = nevals_ul
            evals_LL[fnum, i] = nevals_ll
        end

        println("-----------------------------------------------")

    end

    writedlm("output/accuracy_UL.csv", errors_UL, ',')
    writedlm("output/accuracy_LL.csv", errors_LL, ',')

    writedlm("output/evals_UL.csv", evals_UL, ',')
    writedlm("output/evals_LL.csv", evals_LL, ',')

    ##
    println("Upper level")
    statsToLatex("output/accuracy_UL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")

    println("Lower level")
    statsToLatex("output/accuracy_LL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")


    println("Evaluations")
    statsToLatex("output/evals_UL.csv"; mapping= x->abs.(x))
    println("-------------------------------------------------")


end


test()