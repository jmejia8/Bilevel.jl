# Bilevel Optimization Methods

**Outdated repository.**

[BilevelHeuristics.jl](https://github.com/jmejia8/BilevelHeuristics.jl) is a refactored version of this package.

<hr>

Some methods for Bilevel Optimization for Julia 1.x will be shared here.

## Installation


### Julia 1.0 or Later

First install the **dependencies**: Open the Julia REPL and press `]` to open the Pkg prompt. To add this package, use the add command:

Metaheuristic package:
```
pkg> add https://github.com/jmejia8/Metaheuristics.jl.git
```


BiApprox package:
```
pkg> add https://github.com/jmejia8/BiApprox.git
```



Now, you are able to install Bilevel.jl:
```
pkg> add https://github.com/jmejia8/Bilevel.jl.git
```

## Algorithms

Metaheuristics for Bilevel Optimization.

### SABO: Surrogate-assisted Bilevel Optimization

SABO is a metaheuristic method assisted by a kernel interpolation numerical technique to approximate optimal solutions of a BOP. Two surrogate methods approximate  upper and lower level objective functions on solutions obtained by a population-based algorithm adapted to save upper level objective function evaluations. Some theoretical properties about kernel interpolation are used to study global convergence in some particular BOPs. The empirical results of this approach are analyzed when representative test functions for bilevel optimization are solved. The overall performance provided by the proposal is competitive.

#### Example

```julia
using Bilevel

# upper level objective function
F(x, y) = sum((x + y) .^ 2)

# lower level objective function
f(x, y) = sum((sin.(4π * x) - y) .^ 2)

D = 3 # UL and LL dimension
bounds_ul  = Array([-1 * ones(D) 1 * ones(D)]') # [-1, 1]^D
bounds_ll  = Array([-1 * ones(D) 1 * ones(D)]')

# Information about the bilevel optimization problem
information = Information(F_optimum = 0.0, f_optimum = 0.0)

sabo = SABO(D; information = information)

r = optimize(
    F,
    f,
    bounds_ul,
    bounds_ll,
    sabo
)

display(r)

# optimum (approx)
x, y = r.best_sol.x, r.best_sol.y

```

### QBCA for Bilevel Optimization

QBCA is a population-based metaheuristic approach using Tykhonov regularization and a quasi-Newton method, called Quasi-Newton Bilevel Centers Algorithm (QBCA), to deal with bilevel optimization problems. Tykhonov regularization for bilevel optimization is adopted to handle problems with nonunique lower level solutions. Besides, a quasi-Newton method is adapted to deal with unfeasible solutions in the lower level. The performance of this proposal is assessed by using representative test functions for bilevel optimization.

#### Example

Define the upper level (UL) and lower level (LL) objective functions and run QBCA algorithm.

```julia
using Bilevel

# UL objective function value
F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))

# LL objective function value
f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

bounds = Matrix([-1.0  1.0]')

method = QBCA(F, f, bounds,bounds, s_min=-1, N=20, F_calls_limit = 10000)

result = optimize(method)

best = result.best_sol
```

The `best` structure contains the best solution found. The members of `best` can be accessed as follows:

```julia
println(best.x) # UL vector
println(best.y) # LL vector
println(best.F) # UL objective function value F(x, y)
println(best.f) # LL objective function value f(x, y)
```

### BCA Framework

This framework uses the center of mass concept to generate new solutions at upper
level.


#### Example

Solve the following bilevel optimization problem:

```julia
F(x,y)  = sum((x + y).^2) # upper level
f(x,y)  = sum((x - y).^2) # lower level

# bounds
ul_bounds = [-1 -1; 1 1.0] # -1 <= x[i] <= 1
ll_bounds = [-1 -1; 1 1.0] # -1 <= y[j] <= 1
```

First, specify the options for BCA:

```julia
options = Options(F_calls_limit=1000,
                  f_calls_limit=Int(1e6),
                  debug=true)
```

After that, define the BCA parameters:

```julia

BCA = BCAOperators.BCAFW(N =30, K=7, η_max=2.0) # BCA Framework definition
```

Now, we are able to define the algorithm structure:

```julia

method = Algorithm(BCA;
            initialize! = BCAOperators.initialize!,
            update_state! = BCAOperators.update_state!,
            options = options)
```

Finally, optimize and get results:

```julia
results = optimize(F, f, ul_bounds, ll_bounds, method)

display(results)

```

## TODO
* Documentation
* Bilevel Centers Algorithm
