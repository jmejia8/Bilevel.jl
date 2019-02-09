using Bilevel
using Test

function test1()
    F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
    f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

    bounds = Matrix([-1.0  1.0]')

    method = QBCA(F, f, bounds,bounds, s_min=-1, N=20, F_calls_limit = 10000)

    result = optimize(method)

    best = result.best_sol

    ≈(best.F, 0.8, atol= 1e-1) && ≈(best.f, 0.0, atol=1e-5) 
end

@test test1()