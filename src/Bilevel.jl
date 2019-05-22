module Bilevel

export Algorithm, Information, State, Options, xf_indiv, Problem
export optimize, QBCA

include("structures.jl")
include("externals.jl")
include("display.jl")
include("stop.jl")
include("tools.jl")
include("operators.jl")


include("BCA.jl")
include("QBCA.jl")
include("algorithm.jl")

include("algorithms/Template/procedure.jl")




end # module
