using Bilevel
using Test

F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

function test1()

    bounds = Matrix([-1.0  1.0]')

    method = QBCA(size(bounds, 2); s_min=-1, N=20, F_calls_limit = 10000)

    result = optimize(F, f, bounds,bounds, method)

    best = result.best_sol

    ≈(best.F, 0.8, atol= 1e-1) && ≈(best.f, 0.0, atol=1e-5) 
end

function test2()
    bounds = Matrix([-1.0  1.0]')

    method = QBCA(size(bounds, 2); s_min=-1, N=20, F_calls_limit = 10000, ll_iterations = 10, options = Bilevel.Options(store_convergence = true))

    result = optimize(F, f, bounds,bounds, method)

    best = result.best_sol
    ≈(best.F, 0.8, atol= 1e-1) && ≈(best.f, 0.0, atol=1e-5) && length(result.convergence) > 0
end

function test3()
    method = Algorithm()
    r = optimize((x,y) -> sum(x + y),
             (x,y) -> sum(x + y),
             [0 0; 1 1],
             [0 0; 1 1],
             method
        )

    true
end

@test test1()
@test test2()
@test test3()
