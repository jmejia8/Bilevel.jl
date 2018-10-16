#####################################################
#
#         POULATION METHODS
#
#####################################################

function init_population(N::Int, D::Int, a::Vector{Float64}, b::Vector{Float64})
    # a, b should be D × 1
    # a: lower bounds
    # b: upper bounds
    a'  .+ (b - a)' .* rand(N, D)
end

function init_population(F::Function, f::Function, N::Int, bounds_ul::Matrix, bounds_ll::Matrix)
    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2) 
   
    a_ul, b_ul = bounds_ul[1,:], bounds_ul[2,:]
    a_ll, b_ll = bounds_ll[1,:], bounds_ll[2,:]
   
    X = init_population(N, D_ul, a_ul, b_ul)
    Y = init_population(N, D_ll, a_ll, b_ll)
    
    # infers datatype
    x = X[1,:]
    y = Y[1,:]

    child = generateChild(x, y, F(x, y), f(x, y))
    individual = typeof(child)

    # population array
    population = Array{individual, 1}([])

    # first individual
    push!(population, child)

    for i in 2:N
        x = X[i,:]
        y = Y[i,:]

        child = generateChild(x, y, F(x, y), f(x, y))
        push!(population, child)
    end

    return population
end

#####################################################
#
#         GENRATE NEW SOLUTIONS
#
#####################################################

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Float64, fResult::Float64)
    return xf_indiv(x, y, FResult, fResult)
end

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Tuple{Float64,Array{Float64,1}}, fResult::Tuple{Float64,Array{Float64,1}})
    F, G = FResult
    f, g = fResult
    return xfg_indiv(x, y, F, f, G, g)
end

function generateChild(x::Vector{Float64}, y::Vector{Float64}, FResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}}, fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    F, G, H = FResult
    f, g, h = fResult
    return xfgh_indiv(x, y, F, f, G, g, H, h)
end

function inferType(fVal::Tuple{Float64})
    return xf_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1}})
    return xfg_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    return xfgh_indiv
end
