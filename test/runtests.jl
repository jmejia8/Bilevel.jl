using Bilevel
using Plots
gr(legend=false)

function main()
    F(x, y) = sum(x.^2 + 0.1cos.(4π*x) + y.^2 + 0.1sin.(4π*y))
    f(x, y) = sum((x.^2 + y.^2 .- 1.0).^2)

    bounds = Matrix([-1.0  1.0]')

    println("Optimizing...")
    P, b = optimize(F, f, bounds_ul = bounds, bounds_ll=bounds, η_max=2.0)

    # println("Generating plot...")
    # @gif for i = 1:length(P)
    #     plot(title="$i", xlims=[-1, 1], ylims=[-1, 1])

    #     for indiv = P[i]
    #         scatter!(indiv.x, indiv.y)
    #     end
    # end

end

main()