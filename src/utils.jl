function selectbyproperty(pkproperties, pmin, pmax)
    keep = BitVector(ones(length(pkproperties)))
    if pmin !== nothing
        keep .&= (pmin .<= pkproperties)
    end
    if pmax !== nothing
        keep .&= (pkproperties .<= pmax)
    end
    keep
end

function selectbypeakdistance(pkindices, priority, distance)
    npkindices = length(pkindices)
#    distance_ = ceil.(Int, distance)
    keep = BitVector(ones(npkindices))

    prioritytoposition = fsortperm(priority)
    for i in npkindices:-1:1
        j = prioritytoposition[i]
        (keep[j] == 0) && continue

        k = j-1
        while (1 <= k) && ((pkindices[j]-pkindices[k]) < distance)
            keep[k] = 0
            k -= 1
        end

        k = j+1
        while (k<=npkindices) && ((pkindices[k]-pkindices[j]) < distance)
            keep[k] = 0
            k += 1
        end
    end
    keep
end

function argwlenasexpected(value)
    if value === nothing
        value = -1
    elseif 1 < value
        value = ceil(Int, value)
    else
        throw(ArgumentError("`wlen` must be larger than 1, was $(value)"))
    end
    value
end
