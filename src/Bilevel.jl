module Bilevel

export Algorithm, Information, State, Options, xf_indiv, Problem
export BCAOperators
export optimize, QBCA, SABO

include("types.jl")
include("structures.jl")
include("externals.jl")
include("display.jl")
include("stop.jl")
include("tools.jl")
include("operators.jl")


include("BCA.jl")
include("QBCA.jl")
include("algorithm.jl")

# include("algorithms/Template/Template.jl")
include("algorithms/BCAFW/BCAFW.jl")
include("SABO/SABO.jl")




end # module
