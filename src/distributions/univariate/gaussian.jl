############################################
# GaussianDistribution
############################################
# Description:
#   Encodes a univariate Gaussian distribution.
#   Define (mean (m) or weighted mean (xi))
#   and (variance (V) or precision (W)).
#   These result in the same distribution:
#    dist = GaussianDistribution(m=1.0, V=2.0)
#    dist = GaussianDistribution(m=1.0, W=0.5)
#    dist = GaussianDistribution(xi=0.5, W=0.5)
#    dist = GaussianDistribution(xi=0.5, V=2.0)
############################################

export
    GaussianDistribution,
    ensureMVParametrization!,
    ensureMWParametrization!,
    ensureXiVParametrization!,
    ensureXiWParametrization!,
    ensureParameters!,
    isWellDefined,
    isConsistent

type GaussianDistribution <: UnivariateProbabilityDistribution
    m::Float64   # Mean
    V::Float64   # Variance
    W::Float64   # Precision / weight
    xi::Float64  # Weighted mean: xi=W*m

    function GaussianDistribution(m, V, W, xi)
        self = new(m, V, W, xi)
        isWellDefined(self) || error("Cannot create GaussianDistribution: distribution is underdetermined.")
        if !isnan(V)
            (realmin(Float64) < abs(V) < realmax(Float64)) || error("Cannot create GaussianDistribution: V cannot be 0 or Inf.")
        end
        if !isnan(W)
            (realmin(Float64) < abs(W) < realmax(Float64)) || error("Cannot create GaussianDistribution: W cannot be 0 or Inf.")
        end

        return self
    end
end

function GaussianDistribution(; m::Float64=NaN,
                                V::Float64=NaN,
                                W::Float64=NaN,
                                xi::Float64=NaN)
    return GaussianDistribution(m, V, W, xi)
end

GaussianDistribution() = GaussianDistribution(m=0.0, V=1.0)

vague(::Type{GaussianDistribution}) = GaussianDistribution(m=0.0, V=huge)

function format(dist::GaussianDistribution)
    if !isnan(dist.m) && !isnan(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif !isnan(dist.m) && !isnan(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif !isnan(dist.xi) && !isnan(dist.W)
        return "N(ξ=$(format(dist.xi)), W=$(format(dist.W)))"
    elseif !isnan(dist.xi) && !isnan(dist.V)
        return "N(ξ=$(format(dist.xi)), V=$(format(dist.V)))"
    else
        return "N(underdetermined)"
    end
end

show(io::IO, dist::GaussianDistribution) = println(io, format(dist))

Base.mean(dist::GaussianDistribution) = isProper(dist) ? ensureMDefined!(dist).m : NaN

Base.var(dist::GaussianDistribution) = isProper(dist) ? ensureVDefined!(dist).V : NaN

function isProper(dist::GaussianDistribution)
    if isWellDefined(dist)
        param = isnan(dist.W) ? dist.V : dist.W
        return (realmin(Float64) < param < realmax(Float64))
    else
        return false
    end
end

function sample(dist::GaussianDistribution)
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureMVParametrization!(dist)
    return sqrt(dist.V)*randn() + dist.m
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::GaussianDistribution)
    # Check if dist is not underdetermined
    return ((!isnan(dist.m) || !isnan(dist.xi)) &&
            (!isnan(dist.V) || !isnan(dist.W)))
end

function isConsistent(dist::GaussianDistribution)
    # Check if dist is consistent in case it is overdetermined
    if !isnan(dist.V) && !isnan(dist.W)
        V_W_consistent = false
        if !isApproxEqual(1/dist.V, dist.W)
            return false # V and W are not consistent
        end
    end
    if !isnan(dist.m) && !isnan(dist.xi)
        if !isnan(dist.V)
            if !isApproxEqual(dist.V * dist.xi, dist.m)
                return false
            end
        else
            if !isApproxEqual(dist.W * dist.m, dist.xi)
                return false
            end
        end
    end
    return true # all validations passed
end

function ensureMDefined!(dist::GaussianDistribution)
    # Ensure that dist.m is defined, calculate it if needed.
    dist.m = isnan(dist.m) ? ensureVDefined!(dist).V * dist.xi : dist.m
    return dist
end

function ensureXiDefined!(dist::GaussianDistribution)
    # Ensure that dist.xi is defined, calculate it if needed.
    dist.xi = isnan(dist.xi) ? ensureWDefined!(dist).W * dist.m : dist.xi
    return dist
end

function ensureVDefined!(dist::GaussianDistribution)
    # Ensure that dist.V is defined, calculate it if needed.
    dist.V = isnan(dist.V) ? 1/dist.W : dist.V
    return dist
end

function ensureWDefined!(dist::GaussianDistribution)
    # Ensure that dist.W is defined, calculate it if needed.
    dist.W = isnan(dist.W) ? 1/dist.V : dist.W
    return dist
end

ensureMVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureMDefined!(dist))

ensureMWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureMDefined!(dist))

ensureXiVParametrization!(dist::GaussianDistribution) = ensureVDefined!(ensureXiDefined!(dist))

ensureXiWParametrization!(dist::GaussianDistribution) = ensureWDefined!(ensureXiDefined!(dist))

function ensureParameters!(dist::GaussianDistribution, params::Tuple{Symbol})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return d
end

function ensureParameter!(dist::GaussianDistribution, param::Type{Val{:m}})
    dist.m = isnan(dist.m) ? ensureParameter!(dist, :v).V * dist.xi : dist.m
    return dist
end

function ensureParameter!(dist::GaussianDistribution, param::Type{Val{:V}})
    dist.V = isnan(dist.V) ? 1/dist.W : dist.V
    return dist
end

function ensureParameter!(dist::GaussianDistribution, param::Type{Val{:xi}})
    dist.xi = isnan(dist.xi) ? ensureWDefined!(dist).W * dist.m : dist.xi
    return dist
end

function ensureParameter!(dist::GaussianDistribution, param::Type{Val{:W}})
    dist.W = isnan(dist.W) ? 1/dist.V : dist.W
    return dist
end

function ==(x::GaussianDistribution, y::GaussianDistribution)
    if is(x, y)
        return true
    end
    if !isWellDefined(x) || !isWellDefined(y)
        return false
    end
    # Check m or xi
    if !(isnan(x.m) || isnan(y.m))
        isApproxEqual(x.m, y.m) || return false
    elseif !(isnan(x.xi) || isnan(y.xi))
        isApproxEqual(x.xi, y.xi) || return false
    else
        ensureMDefined!(x); ensureMDefined!(y);
        isApproxEqual(x.m, y.m) || return false
    end

    # Check V or W
    if !(isnan(x.V) || isnan(y.V))
        isApproxEqual(x.V, y.V) || return false
    elseif !(isnan(x.W) || isnan(y.W))
        isApproxEqual(x.W, y.W) || return false
    else
        ensureVDefined!(x); ensureVDefined!(y);
        isApproxEqual(x.V, y.V) || return false
    end

    return true
end

# Converts from DeltaDistribution -> GaussianDistribution
# NOTE: this introduces a small error because the variance is set >0
convert(::Type{GaussianDistribution}, delta::DeltaDistribution{Float64}) = GaussianDistribution(m=delta.m, V=tiny)

convert(::Type{Message{GaussianDistribution}}, msg::Message{DeltaDistribution{Float64}}) = Message(GaussianDistribution(m=msg.payload.m, V=tiny))
