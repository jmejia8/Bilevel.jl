using BilevelBenchmark

include("../src/structures.jl")
include("../src/externals.jl")
include("../src/tools.jl")
include("../src/operators.jl")


include("../src/qca.jl")

function getBilevel(fnum)
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

        P, b, iters, nevals_ul, nevals_ll = optimize(F, f, bounds_ul = bounds_ul, bounds_ll=bounds_ll)
    
        @printf("SMD%d \t F = %e \t f = %e \t ev_ul = %d \t ev_ll = %d\n", fnum, b.F, b.f, nevals_ul, nevals_ll)
    end



end


test()