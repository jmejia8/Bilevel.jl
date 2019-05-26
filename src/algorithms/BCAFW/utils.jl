function fitnessToMass(U)
    fitness = zeros(length(U))
    for i = 1:length(U)
        fitness[i] = U[i].F
    end

    m = minimum(fitness)
    
    m < 0 && (fitness = 2abs(m) .+ fitness)  

    fitness = 2.0maximum(fitness) .- fitness

    return fitness
end

function center!(c, U)
    m = fitnessToMass(U)
    for i = 1:length(U)
        c += m[i] * U[i].x
    end

    return c ./ sum(m)
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

        if a[i] <= c[i] <= b[i]
            x[i] = c[i]
        elseif x[i] < a[i]
            x[i] = a[i]
        else
            x[i] = b[i]
        end

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
