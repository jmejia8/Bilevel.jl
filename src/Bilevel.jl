module Bilevel

export Algorithm, Information, State, Options, xf_indiv, Problem
export BCAOperators, generateChild, xFgh_indiv
export optimize, QBCA, SABO

include("core/types.jl")
include("core/structures.jl")
include("externals.jl")
include("stop.jl")
include("core/tools.jl")
include("operators.jl")


include("BCA.jl")
include("optimize.jl")

# include("algorithms/Template/Template.jl")
include("algorithms/BCAFW/BCAFW.jl")
include("algorithms/SABO/SABO.jl")
include("algorithms/QBCA/QBCA.jl")
include("algorithms/EMBO_MF/EMBO_MF.jl")


include("display.jl")


end # module
