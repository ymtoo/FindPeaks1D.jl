module FindPeaks1D

using SortingLab

export findpeaks1d, peakprominences1d, peakwidths1d

include("utils.jl")

"""
Finds all local maxima in a 1-D signal. The first and the last sample can't be maxima.
"""
function localmaxima1d(x::AbstractVector{T}) where {T<:Real}
    midpts = Vector{Int}(undef, 0)
    leftedges = Vector{Int}(undef, 0)
    rightedges = Vector{Int}(undef, 0)

    i = 2
    imax = length(x)
    while i < imax
        if x[i-1] < x[i]
            iahead = i+1

            while (iahead < imax) && (x[iahead] == x[i])
                iahead += 1
            end

            if x[iahead] < x[i]
                push!(leftedges, i)
                push!(rightedges, iahead-1)
                push!(midpts, (i+iahead-1)÷2)
                i = iahead
            end
        end
        i += 1
    end
    midpts, leftedges, rightedges
end

"""
    findpeaks1d(x; height=nothing, distance=nothing, prominence=nothing, width=nothing, wlen=nothing, relheight=0.5)

Find all local maxima in a 1-D signal with specified `height`, `distance`, `prominence`, `width`.

# Arguments
- `x`: 1-D signal
- `height`: the first element is the minimal and the second, if supplied, is the maximal peak height
- `distance`: the minimal peak distance
- `prominence`: the first element is the minimal and the second, if supplied, is the maximal peak prominence
- `width`: the first element is the minimal and the second, if supplied, is the maximal peak width
- `wlen`: used for calculation of the peak prominence
- `relheight`: used for calculation of peak width

# Returns
Peak indices amd properties

# Examples
```julia-repl
julia> x = [13, 12, 14, 18, 19, 19, 19, 15, 11, 6, 4, 10, 8, 13, 8, 11, 3, 18, 7, 4];

julia> pkindices, pkproperties = findpeaks1d(x)
([6, 12, 14, 16, 18], Dict{String,Any}())

julia> pkindices, pkproperties = findpeaks1d(x, height=11)
([6, 14, 16, 18], Dict{String,Any}("peak_heights" => [19, 13, 11, 18]))

julia> pkindices, pkproperties = findpeaks1d(x, height=11, distance=3)
([6, 14, 18], Dict{String,Any}("peak_heights" => [19, 13, 18]))
```
"""
function findpeaks1d(x::AbstractVector{T}; 
                     height::Union{Nothing,T,NTuple{2,T}}=nothing, 
                     distance::Union{Nothing,I}=nothing, 
                     prominence::Union{Nothing,PT,NTuple{2,PT}}=nothing, 
                     width::Union{Nothing,PT,NTuple{2,PT}}=nothing, 
                     wlen::Union{Nothing,I}=nothing, 
                     relheight::PT=0.5) where {T<:Real,I<:Integer,PT<:AbstractFloat}
    pkindices, leftedges, rightedges = localmaxima1d(x)
    properties = Dict{String,Any}()
    isempty(pkindices) && (return pkindices, properties)

    if height !== nothing
        pkheights = x[pkindices]
        hmin, hmax = height isa Number ? (height, nothing) : height
        keepheight = selectbyproperty(pkheights, hmin, hmax)
        pkindices = pkindices[keepheight]
        properties["peak_heights"] = pkheights
        properties = Dict{String,Any}(key => array[keepheight] for (key, array) in properties)
    end

    if distance !== nothing
        keepdist = selectbypeakdistance(pkindices, x[pkindices], distance)
        pkindices = pkindices[keepdist]
        properties = Dict{String,Any}(key => array[keepdist] for (key, array) in properties)
    end

    if (prominence !== nothing) || (width !== nothing)
        prominences, leftbases, rightbases = peakprominences1d(x, pkindices, wlen)
        properties["prominences"] = prominences
        properties["leftbases"] = leftbases
        properties["rightbases"] = rightbases
    end

    if prominence !== nothing
        pmin, pmax = prominence isa Number ? (prominence, nothing) : prominence
        keepprom = selectbyproperty(prominences, pmin, pmax)
        pkindices = pkindices[keepprom]
        properties = Dict{String,Any}(key => array[keepprom] for (key, array) in properties)
    end

    if width !== nothing
        widths, widthheights, leftips, rightips = peakwidths1d(x, pkindices, relheight, prominences, leftbases, rightbases)
        properties["widths"] = widths
        properties["widthheights"] = widthheights
        properties["leftips"] = leftips
        properties["rightips"] = rightips
        wmin, wmax = width isa Number ? (width, nothing) : width
        keepwidth = selectbyproperty(widths, wmin, wmax)
        pkindices = pkindices[keepwidth]
        properties = Dict{String,Any}(key => array[keepwidth] for (key, array) in properties)
    end

    pkindices, properties
end

"""
    peakprominences1(x, pkindices, wlen=nothing)

Calculate the prominence of each peak in a 1-D signal.

# Arguments
- `x`: 1-D signal
- `pkindices`: peak indices
- `wlen`: a window length in samples to restrict peak finding to a window around the current peak

# Returns
- `prominences`: prominences for each peaks
- `leftbases`, `rightbases`: indices of the left and right of each peak for peaks' bases
"""
function peakprominences1d(x::AbstractVector{T}, 
                           pkindices::AbstractVector{I}, 
                           wlen::Union{Nothing,I}=nothing) where {T<:Real,I<:Integer}
    wlen = argwlenasexpected(wlen)

    prominences = Vector{T}(undef, length(pkindices))
    leftbases = Vector{Int}(undef, length(pkindices))
    rightbases = Vector{Int}(undef, length(pkindices))

    for pknr = 1:length(pkindices)
        pkindex = pkindices[pknr]
        imin = 1
        imax = length(x)
        !(imin <= pkindex <= imax) && throw(ArgumentError("peak $(pkindex) is not a valid index for `x`"))

        if 2 <= wlen
            imin = max(pkindex-wlen÷2, imin)
            imax = min(pkindex+wlen÷2, imax)
        end

        i = leftbases[pknr] = pkindex
        leftmin = x[pkindex]
        while (imin <= i) && (x[i] <= x[pkindex])
            if x[i] < leftmin
                leftmin = x[i]
                leftbases[pknr] = i
            end
            i -= 1
        end

        i = rightbases[pknr] = pkindex
        rightmin = x[pkindex]
        while (i <= imax) && (x[i] <= x[pkindex])
            if x[i] < rightmin
                rightmin = x[i]
                rightbases[pknr] = i
            end
            i += 1
        end

        prominences[pknr] = x[pkindex]-max(leftmin, rightmin)
    end
    prominences, leftbases, rightbases
end

"""
    peakwidths1d(x, pkindices, relheight=0.5, prominencedata=nothing, wlen=nothing)

Calculate the width of each peak in a 1-D signal.

# Arguments
- `x`: 1-D signal
- `pkindices`: peak indices
- `relheight`: relative height with respect to the peak heights and prominences
- `wlen`: a window length in samples to restrict peak finding to a window around the current peak

# Returns
- `widths`: width for each peak in samples
- `widthheights`: height at which the `widths` are evaluated
- `leftips`, `rightips`: interpolated left and right intersection points of a horizontal line at the 
respective `widthheights`
"""
function peakwidths1d(x::AbstractVector{T}, 
                      pkindices::AbstractVector{I}, 
                      relheight::PT=0.5, 
                      prominencedata::Union{Nothing,Tuple}=nothing, 
                      wlen::Union{Nothing,I}=nothing) where {T<:Real,I<:Integer,PT<:AbstractFloat}
    if prominencedata === nothing
        prominencedata = peakprominences1d(x, pkindices, wlen)
    end
    prominences, leftbases, rightbases = prominencedata
    peakwidths1d(x, pkindices, relheight, prominences, leftbases, rightbases)
end
"""
    peakwidths1d(x, pkindices, relheight, prominences, leftbases, rightbases)

Calculate the width of each peak in a 1-D signal.

# Arguments
- `x`: 1-D signal
- `pkindices`: peak indices
- `relheight`: relative height with respect to the peak heights and prominences
- `prominencedata`: output of `peakprominences1d`
- `leftbases`, `rightbases`: indices of the left and right of each peak for peaks' bases

# Returns
- `widths`: width for each peak in samples
- `widthheights`: height at which the `widths` are evaluated
- `leftips`, `rightips`: interpolated left and right intersection points of a horizontal line at the respective `widthheights`
"""
function peakwidths1d(x::AbstractVector{T}, 
                      pkindices::AbstractVector{I}, 
                      relheight::PT, 
                      prominences, 
                      leftbases, 
                      rightbases) where {T<:Real,I<:Integer,PT<:AbstractFloat}

    npkindices = length(pkindices)
    (relheight < 0) && throw(ArgumentError("`relheight` must be greater pr equal to zero"))
    !(npkindices == length(prominences) == length(leftbases) == length(rightbases)) && throw(ArgumentError("arrays in `prominencedata` must have the same length as `pkindices`"))

    widths = Vector{Float64}(undef, npkindices)
    widhtheights = Vector{Float64}(undef, npkindices)
    leftips = Vector{Float64}(undef, npkindices)
    rightips = Vector{Float64}(undef, npkindices)

    for p in 1:npkindices
        imin = leftbases[p]
        imax = rightbases[p]
        pkindex = pkindices[p]
        !(1 <= imin <= pkindex <= imax <= length(x)) && throw(ArgumentError("prominence data is invalid for peak $(pkindex)"))

        height = widhtheights[p] = x[pkindex] - prominences[p] * relheight

        i = pkindex
        while (imin < i) && (height < x[i])
            i -= 1
        end
        leftip = convert(Float64, i)
        if x[i] < height
            leftip += (height-x[i])/(x[i+1]-x[i])
        end
        i = pkindex
        while (i < imax) && height < x[i]
            i += 1
        end
        rightip = convert(Float64, i)
        if x[i] < height
            rightip -= (height-x[i])/(x[i-1]-x[i])
        end

        widths[p] = rightip - leftip
        leftips[p] = leftip
        rightips[p] = rightip
    end
    widths, widhtheights, leftips, rightips
end

end # module