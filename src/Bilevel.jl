module Bilevel

export Algorithm, Information, State, Options, xf_indiv, Problem
export BCAOperators, generateChild, xFgh_indiv
export optimize, QBCA, SABO

include("types.jl")
include("structures.jl")
include("externals.jl")
include("stop.jl")
include("tools.jl")
include("operators.jl")


include("BCA.jl")
include("algorithm.jl")

# include("algorithms/Template/Template.jl")
include("algorithms/BCAFW/BCAFW.jl")
include("SABO/SABO.jl")
include("QBCA/QBCA.jl")


include("display.jl")


end # module
