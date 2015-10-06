############################################
# GaussianDistribution
############################################
# Description:
#   Encodes a univariate or multivariate Gaussian PDF.
#   Define (mean (m) or weighted mean (xi))
#   and (covariance (V) or precision (W)).
#   These result in the same PDF:
#    dist = GaussianDistribution(m=1.0, V=2.0)
#    dist = GaussianDistribution(m=1.0, W=0.5)
#    dist = GaussianDistribution(xi=0.5, W=0.5)
#    dist = GaussianDistribution(xi=0.5, V=2.0)
#   In case of a multivariate distribution,
#   m/xi should be vectors and V/W should be matrices.
############################################

export
    GaussianDistribution,
    ensureMVParametrization!,
    ensureMWParametrization!,
    ensureXiVParametrization!,
    ensureXiWParametrization!,
    isWellDefined,
    isConsistent,
    isProper,
    sample

type GaussianDistribution <: ProbabilityDistribution
    m::Vector{Float64}   # Mean vector
    V::Matrix{Float64}   # Covariance matrix
    W::Matrix{Float64}   # Weight matrix
    xi::Vector{Float64}  # Weighted mean vector: xi=W*m
    function GaussianDistribution(m, V, W, xi)
        (size(m) == size(xi)) || error("Cannot create GaussianDistribution: m and xi should have the same size")
        (size(V) == size(W)) || error("Cannot create GaussianDistribution: V and W should have the same size")

        self = new(m, V, W, xi)
        isWellDefined(self) || error("Cannot create GaussianDistribution, distribution is underdetermined.")
        if isValid(V)
            (maximum(V) < Inf) || error("Cannot create GaussianDistribution, covariance matrix V cannot contain Inf.")
            all(abs(diag(V)) .>= tiny) || error("Cannot create GaussianDistribution, diagonal of covariance matrix V should be non-zero.")
        end
        if isValid(W)
            (maximum(W) < Inf) || error("Cannot create GaussianDistribution, precision matrix W cannot contain Inf.")
            all(abs(diag(W)) .>= tiny) || error("Cannot create GaussianDistribution, diagonal of precision matrix W should be non-zero.")
        end

        return self
    end
end

function GaussianDistribution(; m::Union(Float64,Vector{Float64})=[NaN],
                                V::Union(Float64,Matrix{Float64})=reshape([NaN], 1, 1),
                                W::Union(Float64,Matrix{Float64})=reshape([NaN], 1, 1),
                                xi::Union(Float64,Vector{Float64})=[NaN])
    if typeof(m) <: Vector
        _m = copy(m)
    elseif typeof(m) <: Number
        _m = [m]
    else
        error("m should be a vector or a number")
    end

    if typeof(xi) <: Vector
        _xi = copy(xi)
    elseif typeof(xi) <: Number
        _xi = [xi]
    else
        error("xi should be a vector or a number")
    end

    # Ensure _m and _xi have the same size
    if size(_m) != size(_xi)
        if !isnan(_m[1]) && isnan(_xi[1])
            _xi = fill!(similar(_m), NaN)
        elseif isnan(_m[1]) && !isnan(_xi[1])
            _m = fill!(similar(_xi), NaN)
        else
            error("m and xi should have the same length")
        end
    end

    if typeof(V) <: Matrix
        _V = copy(V)
    elseif typeof(V) <: Number
        _V = fill!(Array(Float64,1,1), V)
    else
        error("V should be a matrix or a number")
    end

    if typeof(W) <: Matrix
        _W = copy(W)
    elseif typeof(W) <: Number
        _W = fill!(Array(Float64,1,1), W)
    else
        error("W should be a matrix or a number")
    end


    # Ensure _V and _W have the same size
    if size(_V) != size(_W)
        if isValid(_V) && !isValid(_W)
            _W = fill!(similar(_V), NaN)
        elseif !isValid(_V) && isValid(_W)
            _V = fill!(similar(_W), NaN)
        else
            error("V and W should have the same size")
        end
    end

    return GaussianDistribution(_m, _V, _W, _xi)
end
GaussianDistribution() = GaussianDistribution(m=0.0, V=1.0)

vague(::Type{GaussianDistribution}) = GaussianDistribution(m=0.0, V=huge)

function format(dist::GaussianDistribution)
    if isValid(dist.m) && isValid(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif isValid(dist.m) && isValid(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif isValid(dist.xi) && isValid(dist.W)
        return "N(ξ=$(format(dist.xi)), W=$(format(dist.W)))"
    else
        return "N(ξ=$(format(dist.xi)), V=$(format(dist.V)))"
    end
end

show(io::IO, dist::GaussianDistribution) = println(io, format(dist))

function Base.mean(dist::GaussianDistribution)
    if isProper(dist)
        return ensureMDefined!(dist).m
    else
        return fill!(similar(dist.m), NaN)
    end
end

function Base.var(dist::GaussianDistribution)
    if isProper(dist)
        return diag(ensureVDefined!(dist).V, 0)
    else
        return fill!(similar(dist.V), NaN)
    end
end

function isProper(dist::GaussianDistribution)
    if isWellDefined(dist)
        param = isValid(dist.V) ? dist.V : dist.W
        if isRoundedPosDef(param)
            return true
        end
    end

    return false
end

function sample(dist::GaussianDistribution)
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureMVParametrization!(dist)
    return (dist.V^0.5)*randn(length(dist.m)) + dist.m
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::GaussianDistribution)
    # Check if dist is not underdetermined
    if ((!isValid(dist.m) && !isValid(dist.xi)) ||
        (!isValid(dist.V) && !isValid(dist.W)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if isValid(getfield(dist, field))
            if dimensions>0
                if maximum(size(getfield(dist, field)))!=dimensions
                    return false
                end
            else
                dimensions = size(getfield(dist, field), 1)
            end
        end
    end
    return true
end

function isConsistent(dist::GaussianDistribution)
    # Check if dist is consistent in case it is overdetermined
    if isValid(dist.V) && isValid(dist.W)
        V_W_consistent = false
        try
           V_W_consistent = isApproxEqual(inv(dist.V), dist.W)
        catch
            try
                V_W_consistent = isApproxEqual(inv(dist.W), dist.V)
            catch
                error("Cannot check consistency of GaussianDistribution because both V and W are non-invertible.")
            end
        end
        if !V_W_consistent
            return false # V and W are not consistent
        end
    end
    if isValid(dist.m) && isValid(dist.xi)
        if isValid(dist.V)
            if isApproxEqual(dist.V * dist.xi, dist.m) == false
                return false # m and xi are not consistent
            end
        else
            if isApproxEqual(dist.W * dist.m, dist.xi) == false
                return false # m and xi are not consistent
            end
        end
    end
    return true # all validations passed
end

function ensureMDefined!(dist::GaussianDistribution)
    # Ensure that dist.m is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    dist.m = !isValid(dist.m) ? ensureVDefined!(dist).V * dist.xi : dist.m
    return dist
end

function ensureXiDefined!(dist::GaussianDistribution)
    # Ensure that dist.xi is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    dist.xi = !isValid(dist.xi) ? ensureWDefined!(dist).W * dist.m : dist.xi
    return dist
end

function ensureVDefined!(dist::GaussianDistribution)
    # Ensure that dist.V is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.V = !isValid(dist.V) ? inv(dist.W) : dist.V
    catch
        error("Cannot calculate V of GaussianDistribution because W is not invertible.")
    end
    return dist
end

function ensureWDefined!(dist::GaussianDistribution)
    # Ensure that dist.W is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.W = !isValid(dist.W) ? inv(dist.V) : dist.W
    catch
        error("Cannot calculate W of GaussianDistribution because V is not invertible.")
    end
    return dist
end

ensureMVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureMDefined!(dist))
ensureMWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureMDefined!(dist))
ensureXiVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureXiDefined!(dist))
ensureXiWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureXiDefined!(dist))

function ==(x::GaussianDistribution, y::GaussianDistribution)
    if is(x, y) return true end
    eps = tiny
    # Check m or xi
    if isValid(x.m) && isValid(y.m)
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    elseif isValid(x.xi) && isValid(y.xi)
        (length(x.xi)==length(x.xi)) || return false
        (maximum(abs(x.xi-y.xi)) < eps) || return false
    else
        ensureMDefined!(x); ensureMDefined!(y);
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    end

    # Check V or W
    if isValid(x.V) && isValid(y.V)
        (length(x.V)==length(x.V)) || return false
        (maximum(abs(x.V-y.V)) < eps) || return false
    elseif isValid(x.W) && isValid(y.W)
        (length(x.W)==length(x.W)) || return false
        (maximum(abs(x.W-y.W)) < eps) || return false
    else
        ensureVDefined!(x); ensureVDefined!(y);
        (length(x.V)==length(x.V)) || return false
        (maximum(abs(x.V-y.V)) < eps) || return false
    end

    return true
end

# Converts from DeltaDistribution -> GaussianDistribution
# NOTE: this introduces a small error because the variance is set >0
convert(::Type{GaussianDistribution}, delta::DeltaDistribution{Float64}) = GaussianDistribution(m=delta.m, V=tiny*eye(length(delta.m)))
convert(::Type{Message{GaussianDistribution}}, msg::Message{DeltaDistribution{Float64}}) = Message(GaussianDistribution(m=msg.payload.m, V=tiny*eye(length(msg.payload.m))))
