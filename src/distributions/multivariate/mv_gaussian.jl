export
    MvGaussianDistribution,
    ensureParameters!,
    isWellDefined,
    isConsistent

"""
Description:

    Encodes a multivariate Gaussian distribution.

Parameters:

    m (“mean”, real vector), V (“variance”, real matrix),
    W (“precision”, real matrix), xi (“weighted mean”, real vector)

Construction:

    MvGaussianDistribution(m=[1.0,3.0], V=[2.0, 0.0; 0.0, 2.0])

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type MvGaussianDistribution{dims} <: MultivariateProbabilityDistribution
    m::Vector{Float64}   # Mean vector
    V::Matrix{Float64}   # Covariance matrix
    W::Matrix{Float64}   # Weight matrix
    xi::Vector{Float64}  # Weighted mean vector: xi=W*m

    function MvGaussianDistribution(m, V, W, xi)
        (size(m) == size(xi)) || error("Cannot create MvGaussianDistribution: m and xi should have the same size")
        (size(V) == size(W)) || error("Cannot create MvGaussianDistribution: V and W should have the same size")
        (length(m) == size(V,1) == size(V,2)) || error("Cannot create MvGaussianDistribution: inconsistent parameter dimensions")

        if isValid(V)
            (maximum(abs(V)) < realmax(Float64)) || error("Cannot create MvGaussianDistribution, covariance matrix V cannot contain Inf.")
            all(abs(diag(V)) .> realmin(Float64)) || error("Cannot create MvGaussianDistribution, diagonal of covariance matrix V should be non-zero.")
        end
        if isValid(W)
            (maximum(abs(W)) < realmax(Float64)) || error("Cannot create MvGaussianDistribution, precision matrix W cannot contain Inf.")
            all(abs(diag(W)) .> realmin(Float64)) || error("Cannot create MvGaussianDistribution, diagonal of precision matrix W should be non-zero.")
        end

        self = new{length(m)}(m, V, W, xi)
        isWellDefined(self) || error("Cannot create MvGaussianDistribution, distribution is underdetermined.")

        return self
    end
end

function MvGaussianDistribution(; m::Vector{Float64}=[NaN],
                                  V::Matrix{Float64}=reshape([NaN], 1, 1),
                                  W::Matrix{Float64}=reshape([NaN], 1, 1),
                                  xi::Vector{Float64}=[NaN])
    # Ensure _m and _xi have the same size
    _m = copy(m)
    _xi = copy(xi)
    if size(_m) != size(_xi)
        if isValid(_m) && !isValid(_xi)
            _xi = fill!(similar(_m), NaN)
        elseif !isValid(_m) && isValid(_xi)
            _m = fill!(similar(_xi), NaN)
        else
            error("m and xi should have the same length")
        end
    end

    # Ensure _V and _W have the same size
    _V = copy(V)
    _W = copy(W)
    if size(_V) != size(_W)
        if isValid(_V) && !isValid(_W)
            _W = fill!(similar(_V), NaN)
        elseif !isValid(_V) && isValid(_W)
            _V = fill!(similar(_W), NaN)
        else
            error("V and W should have the same size")
        end
    end

    return MvGaussianDistribution{length(_m)}(_m, _V, _W, _xi)
end

MvGaussianDistribution() = MvGaussianDistribution(m=zeros(1), V=ones(1,1))

function vague!{dims}(dist::MvGaussianDistribution{dims})
    dist.m = zeros(dims)
    dist.V = huge*eye(dims)
    invalidate!(dist.W)
    invalidate!(dist.xi)
    return dist
end

vague{dims}(::Type{MvGaussianDistribution{dims}}) = MvGaussianDistribution(m=zeros(dims), V=huge*eye(dims))

function format(dist::MvGaussianDistribution)
    if isValid(dist.m) && isValid(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif isValid(dist.m) && isValid(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif isValid(dist.xi) && isValid(dist.W)
        return "N(ξ=$(format(dist.xi)), W=$(format(dist.W)))"
    elseif isValid(dist.xi) && isValid(dist.V)
        return "N(ξ=$(format(dist.xi)), V=$(format(dist.V)))"
    else
        return "N(underdetermined)"
    end
end

show(io::IO, dist::MvGaussianDistribution) = println(io, format(dist))

function Base.mean(dist::MvGaussianDistribution)
    if isProper(dist)
        return ensureParameter!(dist, Val{:m}).m
    else
        return fill!(similar(dist.m), NaN)
    end
end

function Base.cov(dist::MvGaussianDistribution)
    if isProper(dist)
        return ensureParameter!(dist, Val{:V}).V
    else
        return fill!(similar(dist.V), NaN)
    end
end

function Base.var(dist::MvGaussianDistribution)
    if isProper(dist)
        return diag(ensureParameter!(dist, Val{:V}).V)
    else
        return fill!(similar(dist.m), NaN)
    end
end

function isProper(dist::MvGaussianDistribution)
    if isWellDefined(dist)
        param = isValid(dist.V) ? dist.V : dist.W
        if isRoundedPosDef(param)
            return true
        end
    end

    return false
end

function sample(dist::MvGaussianDistribution)
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :V))
    return (dist.V^0.5)*randn(length(dist.m)) + dist.m
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::MvGaussianDistribution)
    # Check if dist is not underdetermined
    if !((isValid(dist.m) || isValid(dist.xi)) &&
         (isValid(dist.V) || isValid(dist.W)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if isValid(getfield(dist, field))
            if dimensions>0
                if maximum(size(getfield(dist, field))) != dimensions
                    return false
                end
            else
                dimensions = size(getfield(dist, field), 1)
            end
        end
    end
    return true
end

function isConsistent(dist::MvGaussianDistribution)
    # Check if dist is consistent in case it is overdetermined
    if isValid(dist.V) && isValid(dist.W)
        V_W_consistent = false
        try
           V_W_consistent = isApproxEqual(inv(cholfact(dist.V)), dist.W)
        catch
            try
                V_W_consistent = isApproxEqual(inv(cholfact(dist.W)), dist.V)
            catch
                error("Cannot check consistency of MvGaussianDistribution because both V and W are non-invertible.")
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

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
# We assume that the distribution is well-defined, otherwise we would've gotten the message upon creating.

function ensureParameters!(dist::MvGaussianDistribution, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

function ensureParameter!(dist::MvGaussianDistribution, param::Type{Val{:m}})
    if !isValid(dist.m)
        dist.m = ensureParameter!(dist, Val{:V}).V * dist.xi
    end
    return dist
end

function ensureParameter!(dist::MvGaussianDistribution, param::Type{Val{:xi}})
    if !isValid(dist.xi)
        dist.xi = ensureParameter!(dist, Val{:W}).W * dist.m
    end
    return dist
end

function ensureParameter!(dist::MvGaussianDistribution, param::Type{Val{:W}})
    if !isValid(dist.W)
        dist.W = inv(cholfact(dist.V))
    end
    return dist
end

function ensureParameter!(dist::MvGaussianDistribution, param::Type{Val{:V}})
    if !isValid(dist.V)
        dist.V = inv(cholfact(dist.W))
    end
    return dist
end

function ==(x::MvGaussianDistribution, y::MvGaussianDistribution)
    if is(x, y)
        return true
    end
    if !isWellDefined(x) || !isWellDefined(y)
        return false
    end
    # Check m or xi
    if isValid(x.m) && isValid(y.m)
        (length(x.m)==length(y.m)) || return false
        isApproxEqual(x.m,y.m) || return false
    elseif isValid(x.xi) && isValid(y.xi)
        (length(x.xi)==length(y.xi)) || return false
        isApproxEqual(x.xi,y.xi) || return false
    else
        ensureParameter!(x, Val{:m}); ensureParameter!(y, Val{:m});
        (length(x.m)==length(y.m)) || return false
        isApproxEqual(x.m,y.m) || return false
    end

    # Check V or W
    if isValid(x.V) && isValid(y.V)
        (length(x.V)==length(y.V)) || return false
        isApproxEqual(x.V,y.V) || return false
    elseif isValid(x.W) && isValid(y.W)
        (length(x.W)==length(y.W)) || return false
        isApproxEqual(x.W,y.W) || return false
    else
        ensureParameter!(x, Val{:V}); ensureParameter!(y, Val{:V});
        (length(x.V)==length(y.V)) || return false
        isApproxEqual(x.V,y.V) || return false
    end

    return true
end

# Convert DeltaDistribution -> MvGaussianDistribution
# NOTE: this introduces a small error because the variance is set >0
convert(::Type{MvGaussianDistribution}, delta::MvDeltaDistribution{Float64}) = MvGaussianDistribution(m=delta.m, V=tiny*eye(length(delta.m)))
convert{TD<:MvDeltaDistribution{Float64}, TG<:MvGaussianDistribution}(::Type{Message{TG}}, msg::Message{TD}) = Message(MvGaussianDistribution(m=msg.payload.m, V=tiny*eye(length(msg.payload.m))))

# Convert GaussianDistribution -> MvGaussianDistribution
convert(::Type{MvGaussianDistribution}, d::GaussianDistribution) = MvGaussianDistribution(m=[d.m], V=d.V*eye(1), W=d.W*eye(1), xi=[d.xi])

# Convert MvGaussianDistribution -> GaussianDistribution
function convert(::Type{GaussianDistribution}, d::MvGaussianDistribution{1})
    GaussianDistribution(m=d.m[1], V=d.V[1,1], W=d.W[1,1], xi=d.xi[1])
end
