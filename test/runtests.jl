using Bilevel
using Test
import Random: seed!

seed!(1)

F(x, y) = sum(x .^ 2 + 0.1 * cos.(4π * x) + y .^ 2 + 0.1 * sin.(4π * y))
f(x, y) = sum((x .^ 2 + y .^ 2 .- 1.0) .^ 2)

function test1()

    bounds = Matrix([-1.0 1.0]')

    method = QBCA(size(bounds, 2); s_min = -1, N = 20)
    method.options.F_calls_limit = 10000

    result = optimize(F, f, bounds, bounds, method)
    best = result.best_sol

    ≈(best.F, 0.8, atol = 1e-1) && ≈(best.f, 0.0, atol = 1e-5)
end

function test2()
    bounds = Matrix([-1.0 1.0]')

    method = QBCA(
        size(bounds, 2);
        s_min = -1,
        N = 20,
        # ll_iterations = 10
    )

    method.options.F_calls_limit = 10000
    method.options.f_calls_limit = 50000
    method.options.store_convergence = true

    result = optimize(F, f, bounds, bounds, method)

    best = result.best_sol
    ≈(best.F, 0.8, atol = 1e-1) &&
        ≈(best.f, 0.0, atol = 1e-5) &&
        length(result.convergence) > 0
end

function test3()

    options = Options(
        F_calls_limit = 1000,
        f_calls_limit = Int(1e6),
        f_tol = 1e-5,
        F_tol = 1e-5,
        debug = false,
        store_convergence = false,
    )


    information = Information(F_optimum = 0.0, f_optimum = 0.0)

    BCA = BCAOperators.BCAFW(N = 30)

    method = Algorithm(
        BCA;
        initialize! = BCAOperators.initialize!,
        update_state! = BCAOperators.update_state!,
        lower_level_optimizer = BCAOperators.lower_level_optimizer,
        is_better = BCAOperators.is_better,
        stop_criteria = BCAOperators.stop_criteria,
        options = options,
        information = information,
    )

    r = optimize(
        (x, y) -> sum((x + y) .^ 2),
        (x, y) -> sum((sin.(4π * x) - y) .^ 2),
        [-1 -1; 1 1.0],
        [-1 -1; 1 1.0],
        method,
    )



    true
end

function testSABO()
    D = 5
    bounds = Array([-1 * ones(D) 1 * ones(D)]')

    information = Information(F_optimum = 0.0, f_optimum = 0.0)


    # F(x, y) = sum((x + y) .^ 2)
    # f(x, y) = sum((sin.(4π * x) - y) .^ 2)


    F(x, y) = sum( (x[1:end-3] .- 0.2) .^2 ) + sum( (y[end-2:end] - x[end-2:end]) .^2 )
    f(x, y) = (x[1] .- 0.2)^2 + sum( (y[1:end-3] .- 0.2) .^2 ) + sum( (y[end-2:end] - x[end-2:end]) .^2 )

    method = SABO(D, N = 20; information = information)
    r = optimize(
        F, f,
        bounds,
        bounds,
        method
    )
    display(r)
    #r.best_sol.F < 1e-2 && r.best_sol.f < 1e-3
    true

end

function test_mo()
    FF(x, y) = ( [x[1], x[2]], [ sum(x - y), sum(y) ], [prod(x + y), sum(x)] )
    ff(x, y) = ( [x[1]-y[1], x[2] - y[1]], [ sum(x - y) ], [prod(x + y)] )

    x = rand(2)
    y = rand(2)

    sol = generateChild(x, y, FF(x, y), ff(x, y))
    true
end


# @test test1()
# @test test2()
# @test test3()
@test testSABO()
# @test test_mo()
