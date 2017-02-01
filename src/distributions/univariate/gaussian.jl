export
    Gaussian,
    ensureParameters!,
    isWellDefined,
    isConsistent

"""
Description:

    Encodes a univariate Gaussian distribution.
    The Gaussian distribution can be parametrized in multiple ways.
    Depending on the application, a specific parametrization might be attractive from a computational point of view.
    The following combinations are valid: (m,V), (m,W), (xi,V), (xi,W).

Parameters:

    Real scalars m (mean), V (variance), W (precision), xi (weighted mean)

Construction:

    These result in the same distribution:
    Gaussian(m=1.0, V=2.0)
    Gaussian(m=1.0, W=0.5)
    Gaussian(xi=0.5, W=0.5)
    Gaussian(xi=0.5, V=2.0)

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type Gaussian <: Univariate
    m::Float64   # Mean
    V::Float64   # Variance
    W::Float64   # Precision / weight
    xi::Float64  # Weighted mean: xi=W*m

    function Gaussian(m, V, W, xi)
        self = new(m, V, W, xi)
        isWellDefined(self) || error("Cannot create Gaussian: distribution is underdetermined.")
        if !isnan(V)
            (realmin(Float64) < abs(V) < realmax(Float64)) || error("Cannot create Gaussian: V cannot be 0 or Inf.")
        end
        if !isnan(W)
            (realmin(Float64) < abs(W) < realmax(Float64)) || error("Cannot create Gaussian: W cannot be 0 or Inf.")
        end

        return self
    end
end

function Gaussian(; m::Float64=NaN,
                    V::Float64=NaN,
                    W::Float64=NaN,
                    xi::Float64=NaN)
    if isnan(m) && isnan(V) && isnan(W) && isnan(xi)
        m = 0.0
        V = 1.0
    end

    return Gaussian(m, V, W, xi)
end

function pdf(dist::Gaussian, x::Float64)
    ensureParameter!(dist, Val{:m})
    if !isnan(dist.W)
        C = sqrt(dist.W) / sqrt(2.0*pi)
        return C * exp(-0.5 * (x-dist.m)^2 * dist.W)
    elseif !isnan(dist.V)
        C = 1.0 / sqrt(2.0*pi*dist.V)
        return C * exp(-0.5 * (x-dist.m)^2 / dist.V)
    else
        error("Cannot evaluate pdf for underdetermined Gaussian")
    end
end

function logpdf(dist::Gaussian, x::Float64)
    ensureParameter!(dist, Val{:m})
    if !isnan(dist.W)
        C = 0.5*log(dist.W) - 0.5*log(2.0*pi)
        return C - 0.5*dist.W*(x - dist.m)^2
    elseif !isnan(dist.V)
        C = -0.5*log(dist.V) - 0.5*log(2.0*pi)
        return C - 0.5*(1/dist.V)*(x - dist.m)^2
    else
        error("Cannot evaluate logpdf for underdetermined Gaussian")
    end
end

function vague!(dist::Gaussian)
    dist.m = 0.0
    dist.V = huge
    dist.W = NaN
    dist.xi = NaN
    return dist
end

function format(dist::Gaussian)
    if !isnan(dist.m) && !isnan(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif !isnan(dist.m) && !isnan(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif !isnan(dist.xi) && !isnan(dist.W)
        ensureParameters!(dist, (:m,))
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif !isnan(dist.xi) && !isnan(dist.V)
        ensureParameters!(dist, (:m,))
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    else
        return "N(underdetermined)"
    end
end

show(io::IO, dist::Gaussian) = println(io, format(dist))

unsafeMean(dist::Gaussian) = ensureParameter!(dist, Val{:m}).m # unsafe mean

unsafeCov(dist::Gaussian) = ensureParameter!(dist, Val{:V}).V # unsafe covariance

unsafeVar(dist::Gaussian) = unsafeCov(dist) # unsafe variance

function prod!(x::Gaussian, y::Gaussian, z::Gaussian=Gaussian())
    # Multiplication of 2 Gaussian PDFs: p(z) = p(x) * p(y)
    ensureParameters!(x, (:xi, :W))
    ensureParameters!(y, (:xi, :W))
    z.m = NaN
    z.V = NaN
    z.W  = x.W + y.W
    z.xi = x.xi + y.xi

    return z
end

@symmetrical function prod!(x::Gaussian, y::Delta{Float64}, z::Delta{Float64}=Delta())
    # Multiplication of Gaussian PDF with a Delta
    z.m = y.m

    return z
end

@symmetrical function prod!(x::Gaussian, y::Delta{Float64}, z::Gaussian)
    # Multiplication of Gaussian PDF with Delta, force result to be Gaussian
    z.m = y.m
    z.V = tiny
    z.xi = z.W = NaN

    return z
end

@symmetrical function prod!(::Void, y::Delta{Float64}, z::Gaussian)
    # Multiplication of an unknown with Delta, force result to be Gaussian
    z.m = y.m
    z.V = tiny
    z.xi = z.W = NaN

    return z
end

function isProper(dist::Gaussian)
    if isWellDefined(dist)
        param = isnan(dist.W) ? dist.V : dist.W
        return (realmin(Float64) < param < realmax(Float64))
    else
        return false
    end
end

function sample(dist::Gaussian)
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :V))
    return sqrt(dist.V)*randn() + dist.m
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::Gaussian)
    # Check if dist is not underdetermined
    return ((!isnan(dist.m) || !isnan(dist.xi)) &&
            (!isnan(dist.V) || !isnan(dist.W)))
end

function isConsistent(dist::Gaussian)
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

function ensureParameters!(dist::Gaussian, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
# We assume that the distribution is well-defined, otherwise we would've gotten the message upon creating.

function ensureParameter!(dist::Gaussian, param::Type{Val{:m}})
    if isnan(dist.m)
        dist.m = ensureParameter!(dist, Val{:V}).V * dist.xi
    end
    return dist
end

function ensureParameter!(dist::Gaussian, param::Type{Val{:V}})
    if isnan(dist.V)
        dist.V = 1/dist.W
    end
    return dist
end

function ensureParameter!(dist::Gaussian, param::Type{Val{:xi}})
    if isnan(dist.xi)
        dist.xi = ensureParameter!(dist, Val{:W}).W * dist.m
    end
    return dist
end

function ensureParameter!(dist::Gaussian, param::Type{Val{:W}})
    if isnan(dist.W)
        dist.W = 1/dist.V
    end
    return dist
end

function ==(x::Gaussian, y::Gaussian)
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
        ensureParameter!(x, Val{:m}); ensureParameter!(y, Val{:m});
        isApproxEqual(x.m, y.m) || return false
    end

    # Check V or W
    if !(isnan(x.V) || isnan(y.V))
        isApproxEqual(x.V, y.V) || return false
    elseif !(isnan(x.W) || isnan(y.W))
        isApproxEqual(x.W, y.W) || return false
    else
        ensureParameter!(x, Val{:V}); ensureParameter!(y, Val{:V});
        isApproxEqual(x.V, y.V) || return false
    end

    return true
end

# Converts from Delta -> Gaussian
# NOTE: this introduces a small error because the variance is set >0
convert(::Type{Gaussian}, flt::Float64) = Gaussian(m=flt, V=tiny)

convert(::Type{Gaussian}, delta::Delta{Float64}) = Gaussian(m=delta.m, V=tiny)

convert(::Type{Message{Gaussian}}, msg::Message{Delta{Float64}}) = Message(Gaussian(m=msg.payload.m, V=tiny))

# Entropy functional
function differentialEntropy(dist::Gaussian)
    ensureParameters!(dist, (:m, :V))

    return  0.5*log(dist.V) +
            0.5*log(2*pi) +
            0.5
end