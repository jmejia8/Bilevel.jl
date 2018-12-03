≺(a, b) = Selection(b, a)

function Selection(Old::xf_indiv, New::xf_indiv)
    return is_better(New, Old)
end

# Deb rules (selection)
function Selection(Old::xfgh_indiv, New::xfgh_indiv)

    old_vio = violationsSum(Old.g, Old.h) + violationsSum(Old.G, Old.H)
    new_vio = violationsSum(New.g, New.h) + violationsSum(New.G, New.H)

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    return is_better(New, Old)
end

# Deb rules (selection)
function Selection(Old::xfg_indiv, New::xfg_indiv)
    old_vio = violationsSum(Old.g, []) + violationsSum(Old.g, [])
    new_vio = violationsSum(New.g, []) + violationsSum(Old.G, [])

    if new_vio < old_vio 
        return true
    elseif new_vio > old_vio 
        return false
    end

    return is_better(New, Old)
end

function is_better(S1, S2)
    # S1 is better than S2

    return 0.1(S1.F) + (S1.f) < 0.1(S2.F) + (S2.f)
    # return abs(S1.F) + abs(S1.f) < abs(S2.F) + abs(S2.f)
end

function getBest(Population::Array)
    best = Population[1]

    for i = 2:length(Population)
        if Population[i] ≺ best
            best = Population[i]
        end
    end

    return best
end

function getWorstInd(Population::Array)
    worst = 1

    for i = 2:length(Population)
        if Population[worst] ≺ Population[i]
            worst = i
        end
    end

    return worst
end

getWorst(Population::Array) = Population[getWorstInd(Population)]
