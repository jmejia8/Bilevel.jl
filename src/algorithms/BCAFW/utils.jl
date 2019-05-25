function center!(c, U)
    m = 0.0
    for u = U
        c += u.f * u.x
        m += u.f
    end

    return c ./ m
end

function correctSol!(x, a, b)
    # Correct solution
    for i = 1:length(x)
        if a[i] <= x[i] <= b[i]
            continue
        end

        if x[i] < a[i]
            x[i] = a[i]
        else
            x[i] = b[i]
        end
    end
    
    return x
end

function correctSol!(x,c, a, b)
    # Correct solution
    for i = 1:length(x)
        if a[i] <= x[i] <= b[i]
            continue
        end

        x[i] = c[i]
    end
    
    return x
end

function worst(U,is_better)
    u_worst = U[1]
    for u in U
        u_worst = is_better(u_worst, u) ? u : u_worst
    end
    u_worst.x
end

function getWorstInd(population, is_better)
    worst = 1

    for i = 2:length(population)
        if is_better(population[worst],  population[i])
            worst = i
        end
    end

    return worst
end