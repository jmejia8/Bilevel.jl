# Bilevel Optimization Methods

Some methods for Bilevel Optimization for Julia 1.* will be shared here.

## Installation

### Julia 0.7 or Later

Open the Julia REPL and press `]` to open the Pkg prompt. To add this package, use the add command:
```
pkg> add https://github.com/jmejia8/Bilevel.jl.git
```

## Algorithms

Metaheuristics for Bilevel Optimization.

### QBCA for Bilevel Optimization
This work presents a population-based metaheuristic approach using Tykhonov regularization and a quasi-Newton method, called Quasi-Newton Bilevel Centers Algorithm (QBCA), to deal with bilevel optimization problems. Tykhonov regularization for bilevel optimization is adopted to handle problems with nonunique lower level solutions. Besides, a quasi-Newton method is adapted to deal with unfeasible solutions in the lower level. The performance of this proposal is assessed by using representative test functions for bilevel optimization.

#### Example

```julia
using Bilevel
F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

bounds = Matrix([-1.0  1.0]')

method = QBCA(F, f, bounds,bounds, s_min=-1, N=20, F_calls_limit = 10000)

result = optimize(method)

best = result.best_sol
```

## TODO
* Documentation
* Bilevel Centers Algorithm 