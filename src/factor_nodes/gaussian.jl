export Gaussian

abstract Gaussian <: SoftFactor

Univariate(family::Type{Gaussian}; kwargs...) = Univariate{family}(Dict(kwargs))
Multivariate(family::Type{Gaussian}; kwargs...) = Multivariate{family, size(kwargs[1][2])[1]}(Dict(kwargs))

vague(::Type{Univariate{Gaussian}}) = Univariate(Gaussian, m=0.0, v=huge)
vague{dims}(::Type{Multivariate{Gaussian, dims}}) = Multivariate(Gaussian, m=zeros(dims), v=huge*diageye(dims))

unsafeMean(dist::ProbabilityDistribution{Gaussian}) = deepcopy(ensureParameter!(dist, Val{:m}).params[:m]) # unsafe mean

unsafeVar(dist::Univariate{Gaussian}) = ensureParameter!(dist, Val{:v}).params[:v] # unsafe variance
unsafeVar(dist::Multivariate{Gaussian}) = diag(ensureParameter!(dist, Val{:v}).params[:v])

unsafeCov(dist::Univariate{Gaussian}) = ensureParameter!(dist, Val{:v}).params[:v] # unsafe covariance
unsafeCov(dist::Multivariate{Gaussian}) = deepcopy(ensureParameter!(dist, Val{:v}).params[:v])

function isProper(dist::Univariate{Gaussian})
    if isWellDefined(dist)
        if isValid(dist, :w)
            param = dist.params[:w]
        elseif isValid(dist, :v)
            param = dist.params[:v]
        else
            return false
        end

        return (realmin(Float64) < param < realmax(Float64))
    end
    
    return false
end

function isProper(dist::Multivariate{Gaussian})
    if isWellDefined(dist)
        if isValid(dist, :w)
            param = dist.params[:w]
        elseif isValid(dist, :v)
            param = dist.params[:v]
        else
            return false
        end
        
        return isRoundedPosDef(param)
    end
    
    return false
end

function isWellDefined(dist::Univariate{Gaussian})
    # Check if dist is not underdetermined
    location_valid = isValid(dist, :m) || isValid(dist, :xi)
    scale_valid    = isValid(dist, :v) || isValid(dist, :w)

    return location_valid && scale_valid
end

function isWellDefined(dist::Multivariate{Gaussian})
    # Check if dist is not underdetermined
    location_valid = isValid(dist, :m) || isValid(dist, :xi)
    scale_valid    = isValid(dist, :v) || isValid(dist, :w)

    if !location_valid || !scale_valid
        return false
    end

    dimensions=0
    for field in [:m, :xi, :v, :w]
        if isValid(dist, field)
            if dimensions>0
                if maximum(size(dist.params[field])) != dimensions
                    return false
                end
            else
                dimensions = size(dist.params[field], 1)
            end
        end
    end

    return true
end

function sample(dist::Univariate{Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :v))
    return sqrt(dist.params[:v])*randn() + dist.params[:m]
end

function sample{dims}(dist::Multivariate{Gaussian, dims})
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :v))
    return chol(dist.params[:v])' *randn(dims) + dist.params[:m]
end

function prod!( x::Univariate{Gaussian},
                y::Univariate{Gaussian},
                z::Univariate{Gaussian}=Univariate(Gaussian, xi=0.0, w=1.0))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    z.params[:m] = NaN
    z.params[:v] = NaN
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

function prod!{dims}(x::Multivariate{Gaussian, dims},
                     y::Multivariate{Gaussian, dims},
                     z::Multivariate{Gaussian, dims}=Multivariate(Gaussian, xi=zeros(dims), w=diageye(dims)))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    invalidate!(z, :m)
    invalidate!(z, :v)
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

function ensureParameters!(dist::ProbabilityDistribution, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:m}})
    if !isValid(dist, :m)
        dist.params[:m] = ensureParameter!(dist, Val{:v}).params[:v] * dist.params[:xi]
    end
    return dist
end

function ensureParameter!(dist::Univariate{Gaussian}, param::Type{Val{:v}})
    if !isValid(dist, :v)
        dist.params[:v] = 1/dist.params[:w]
    end
    return dist
end

function ensureParameter!(dist::Multivariate{Gaussian}, param::Type{Val{:v}})
    if !isValid(dist, :v)
        dist.params[:v] = cholinv(dist.params[:w])
    end
    return dist
end

function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:xi}})
    if !isValid(dist, :xi)
        dist.params[:xi] = ensureParameter!(dist, Val{:w}).params[:w] * dist.params[:m]
    end
    return dist
end

function ensureParameter!(dist::Univariate{Gaussian}, param::Type{Val{:w}})
    if !isValid(dist, :w)
        dist.params[:w] = 1/dist.params[:v]
    end
    return dist
end

function ensureParameter!(dist::Multivariate{Gaussian}, param::Type{Val{:w}})
    if !isValid(dist, :w)
        dist.params[:w] = cholinv(dist.params[:v])
    end
    return dist
end

function ==(t::Univariate{Gaussian}, u::Univariate{Gaussian})
    if is(t, u)
        return true
    end
    ensureParameters!(t, (:xi, :w))
    ensureParameters!(u, (:xi, :w))
    return isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

function ==(t::Multivariate{Gaussian}, u::Multivariate{Gaussian})
    if is(t, u)
        return true
    end
    ensureParameters!(t, (:xi, :w))
    ensureParameters!(u, (:xi, :w))
    return isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

# Entropy functional
function differentialEntropy(dist::Univariate{Gaussian})
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy{dims}(dist::Multivariate{Gaussian, dims})
    return  0.5*log(det(unsafeCov(dist))) +
            (dims/2)*log(2*pi) +
            (dims/2)
end