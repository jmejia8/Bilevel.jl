using Bilevel
using Test

function test1()
    F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
    f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

    bounds = Matrix([-1.0  1.0]')

    method = QBCA(F, f, bounds,bounds)

    result = optimize(method)

    best = result.best_sol
    println(best.F, " ", best.f)

    (best.F - 0.82) < 01e-2 && best.f ≈ 0.0 
end

@test test1()