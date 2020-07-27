#####################################################
#
#         POULATION METHODS
#
#####################################################

function init_population(N::Int, D::Int, a::Vector{Float64}, b::Vector{Float64})
    # a, b should be D × 1
    # a: lower bounds
    # b: upper bounds
    a' .+ (b - a)' .* rand(N, D)
end

function init_population(
    F::Function,
    f::Function,
    N::Int,
    bounds_ul::Matrix,
    bounds_ll::Matrix,
)
    D_ul, D_ll = size(bounds_ul, 2), size(bounds_ll, 2)

    a_ul, b_ul = bounds_ul[1, :], bounds_ul[2, :]
    a_ll, b_ll = bounds_ll[1, :], bounds_ll[2, :]

    X = init_population(N, D_ul, a_ul, b_ul)
    Y = init_population(1, D_ll, a_ll, b_ll)

    # infers datatype
    x = X[1, :]
    y = Y[1, :]

    nevals_ll = 0

    child = generateChild(x, y, F(x, y), f(x, y))
    individual = typeof(child)

    nevals_ll += 1

    # population array
    population = Array{individual,1}([])

    # first individual
    push!(population, child)

    for i = 2:N
        x = X[i, :]

        opt = Metaheuristics.Options(f_calls_limit = 1000D_ll)
        eca = Metaheuristics.ECA(;
            K = 3,
            N = 3D_ll,
            p_bin = 0.0,
            p_exploit = 2.0,
            options = opt,
        )
        res = Metaheuristics.optimize(z -> f(x, z), bounds_ll, eca)

        # y, fy = Metaheuristics.eca( z-> f(x, z), D_ll;
        #                                 limits=bounds_ll,
        #                                 K = 3, N = 3D_ll,
        #                                 showResults=false,
        #                                 p_bin = 0,
        #                                 p_exploit = 2,
        #                                 canResizePop=false,
        #                                 max_evals=1000D_ll)
        y = res.best_sol.x
        fy = res.best_sol.f

        y = correct(y, bounds_ll)

        r = Optim.optimize(
            z -> f(x, z),
            bounds_ll[1, :],
            bounds_ll[2, :],
            y,
            Optim.Fminbox(Optim.BFGS()),
        )
        y = r.minimizer

        nevals_ll += res.f_calls + r.f_calls + 1

        child = generateChild(x, y, F(x, y), f(x, y))
        push!(population, child)
    end

    return population, nevals_ll
end

#####################################################
#
#         GENRATE NEW SOLUTIONS
#
#####################################################

function generateChild(
    x::Vector{Float64},
    y,
    FResult::Float64,
    fResult::Float64,
)
    return xf_indiv(x, y, FResult, fResult)
end

function generateChild(
    x::Vector{Float64},
    y,
    FResult::Tuple{Float64,Array{Float64,1}},
    fResult::Tuple{Float64,Array{Float64,1}},
)
    F, G = FResult
    f, g = fResult
    return xfg_indiv(x, y, F, f, G, g)
end

function generateChild(
    x::Vector{Float64},
    y,
    FResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}},
    fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}},
)
    F, G, H = FResult
    f, g, h = fResult
    return xfgh_indiv(x, y, F, f, G, g, H, h)
end

function generateChild(
    x::Vector{Float64},
    y,
    FResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}},
    fResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}},
)
    F, G, H = FResult
    f, g, h = fResult

    return xFgh_indiv(x, y, F, f, G, g, H, h)
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

function inferType(fVal::Tuple{Array{Float64},Array{Float64,1},Array{Float64,1}})
    return xFgh_indiv
end

############################################
function correct(x::Vector{Float64}, bounds::Array)
    # Correct solution

    for i = 1:length(x)
        if !(bounds[1, i] <= x[i] <= bounds[2, i])
            # print(x[i], "  --->   ")
            x[i] = bounds[1, i] + (bounds[2, i] - bounds[1, i]) * rand()
            # println(x[i])
        end
    end

    return x
end
