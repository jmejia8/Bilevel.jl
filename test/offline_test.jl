using BilevelBenchmark

include("../src/structures.jl")
include("../src/externals.jl")
include("../src/tools.jl")
include("../src/operators.jl")


include("../src/qca.jl")

function test()
    fnum = 3
    D_ul = D_ll = 5

    bounds_ul, bounds_ll = bilevel_ranges(D_ul, D_ll, fnum)

    # leader
    F(x::Array{Float64}, y::Array{Float64}) =  bilevel_leader(x, y, fnum)

    # follower
    f(x::Array{Float64}, y::Array{Float64}) = bilevel_follower(x, y, fnum)


    println("Optimizing...")
    P, b = optimize(F, f, bounds_ul = bounds_ul, bounds_ll=bounds_ll)


    println("f: ", b.f)
    println("F: ", b.F)

end


test()