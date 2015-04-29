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
    sample

type GaussianDistribution <: ProbabilityDistribution
    m::Vector{Float64}   # Mean vector
    V::Matrix{Float64}   # Covariance matrix
    W::Matrix{Float64}   # Weight matrix
    xi::Vector{Float64}  # Weighted mean vector: xi=W*m
    function GaussianDistribution(m, V, W, xi)
        self = new(m, V, W, xi)
        isWellDefined(self) || error("Cannot create GaussianDistribution, distribution is underdetermined.")
        if valid(V)
            (maximum(V) < Inf) || error("Cannot create GaussianDistribution, covariance matrix V cannot contain Inf.")
            (isRoundedPosDef(V)) || error("Cannot create GaussianDistribution, covariance matrix V should be positive definite.")
        end
        if valid(W)
            (maximum(W) < Inf) || error("Cannot create GaussianDistribution, precision matrix W cannot contain Inf.")
            (isRoundedPosDef(W)) || error("Cannot create GaussianDistribution, precision matrix W should be positive definite.")
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
        _m = [NaN]
    end

    if typeof(xi) <: Vector
        _xi = copy(xi)
    elseif typeof(xi) <: Number
        _xi = [xi]
    else
        _xi = [NaN]
    end

    if typeof(V) <: Matrix
        _V = copy(V)
    elseif typeof(V) <: Number
        _V = fill!(Array(Float64,1,1), V)
    else
        _V = reshape([NaN], 1, 1)
    end

    if typeof(W) <: Matrix
        _W = copy(W)
    elseif typeof(W) <: Number
        _W = fill!(Array(Float64,1,1), W)
    else
        _W = reshape([NaN], 1, 1)
    end

    return GaussianDistribution(_m, _V, _W, _xi)
end
GaussianDistribution() = GaussianDistribution(m=0.0, V=1.0)

vague(::Type{GaussianDistribution}) = GaussianDistribution(m=0.0, V=huge())

function format(dist::GaussianDistribution)
    if valid(dist.m) && valid(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif valid(dist.m) && valid(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif valid(dist.xi) && valid(dist.W)
        return "N(ξ=$(format(dist.xi)), W=$(format(dist.W)))"
    else
        return "N(ξ=$(format(dist.xi)), V=$(format(dist.V)))"
    end        
end

show(io::IO, dist::GaussianDistribution) = println(io, format(dist))

Base.mean(dist::GaussianDistribution) = ensureMDefined!(dist).m
Base.var(dist::GaussianDistribution) = diag(ensureVDefined!(dist).V, 0)

function sample(dist::GaussianDistribution)
    ensureMVParametrization!(dist)
    return DeltaDistribution((dist.V^0.5)*randn(length(dist.m)) + dist.m)
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::GaussianDistribution)
    # Check if dist is not underdetermined
    if ((!valid(dist.m) && !valid(dist.xi)) ||
        (!valid(dist.V) && !valid(dist.W)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if valid(getfield(dist, field))
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
    if valid(dist.V) && valid(dist.W)
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
    if valid(dist.m) && valid(dist.xi)
        if valid(dist.V)
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
    dist.m = !valid(dist.m) ? ensureVDefined!(dist).V * dist.xi : dist.m
    return dist
end
function ensureXiDefined!(dist::GaussianDistribution)
    # Ensure that dist.xi is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    dist.xi = !valid(dist.xi) ? ensureWDefined!(dist).W * dist.m : dist.xi
    return dist
end
function ensureVDefined!(dist::GaussianDistribution)
    # Ensure that dist.V is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.V = !valid(dist.V) ? inv(dist.W) : dist.V
    catch
        error("Cannot calculate V of GaussianDistribution because W is not invertible.")
    end
    return dist
end
function ensureWDefined!(dist::GaussianDistribution)
    # Ensure that dist.W is defined, calculate it if needed.
    # An underdetermined dist will throw an exception, we assume dist is well defined.
    try
        dist.W = !valid(dist.W) ? inv(dist.V) : dist.W
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
    eps = tiny()
    # Check m or xi
    if valid(x.m) && valid(y.m)
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    elseif valid(x.xi) && valid(y.xi)
        (length(x.xi)==length(x.xi)) || return false
        (maximum(abs(x.xi-y.xi)) < eps) || return false
    else
        ensureMDefined!(x); ensureMDefined!(y);
        (length(x.m)==length(x.m)) || return false
        (maximum(abs(x.m-y.m)) < eps) || return false
    end

    # Check V or W
    if valid(x.V) && valid(y.V)
        (length(x.V)==length(x.V)) || return false
        (maximum(abs(x.V-y.V)) < eps) || return false
    elseif valid(x.W) && valid(y.W)
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
convert{T<:Real}(::Type{GaussianDistribution}, delta::DeltaDistribution{T}) = GaussianDistribution(m=delta.m, V=tiny())
convert{T<:Real}(::Type{GaussianDistribution}, delta::DeltaDistribution{Vector{T}}) = GaussianDistribution(m=deepcopy(delta.m), V=tiny()*eye(length(delta.m)))
convert{T<:Real}(::Type{Message{GaussianDistribution}}, msg::Message{DeltaDistribution{T}}) = Message(GaussianDistribution(m=msg.payload.m, V=tiny()))
convert{T<:Real}(::Type{Message{GaussianDistribution}}, msg::Message{DeltaDistribution{Vector{T}}}) = Message(GaussianDistribution(m=deepcopy(msg.payload.m), V=tiny()*eye(length(msg.payload.m))))
