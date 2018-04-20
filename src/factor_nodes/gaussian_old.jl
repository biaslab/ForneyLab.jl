
slug(::Type{Gaussian}) = "𝒩"

function format{V<:VariateType}(dist::ProbabilityDistribution{V, Gaussian})
    ensureParameters!(dist, (:m, :v))
    return "$(slug(Gaussian))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"
end

ProbabilityDistribution(::Type{Univariate}, ::Type{Gaussian}; m::Number=NaN, v::Number=NaN, w::Number=NaN, xi::Number=NaN) = ProbabilityDistribution{Univariate, Gaussian}(Dict(:m=>m, :v=>v, :w=>w, :xi=>xi))
ProbabilityDistribution(::Type{Gaussian}; m::Number=NaN, v::Number=NaN, w::Number=NaN, xi::Number=NaN) = ProbabilityDistribution{Univariate, Gaussian}(Dict(:m=>m, :v=>v, :w=>w, :xi=>xi))
ProbabilityDistribution(::Type{Multivariate}, ::Type{Gaussian}; m::Vector=[NaN], v::AbstractMatrix=[NaN].', w::AbstractMatrix=[NaN].', xi::Vector=[NaN]) = ProbabilityDistribution{Multivariate, Gaussian}(Dict(:m=>m, :v=>v, :w=>w, :xi=>xi))

dims(dist::ProbabilityDistribution{Univariate, Gaussian}) = 1
function dims(dist::ProbabilityDistribution{Multivariate, Gaussian})
    if isValid(dist, :m)
        return length(dist.params[:m])
    elseif isValid(dist, :xi)
        return length(dist.params[:xi])
    else
        error("Invalid distribution $(dist)")
    end
end

vague(::Type{Gaussian}) = ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=huge)
vague(::Type{Gaussian}, dims::Int64) = ProbabilityDistribution(Multivariate, Gaussian, m=zeros(dims), v=huge*diageye(dims))

unsafeMean{T<:VariateType}(dist::ProbabilityDistribution{T, Gaussian}) = deepcopy(ensureParameter!(dist, Val{:m}).params[:m]) # unsafe mean

unsafeVar(dist::ProbabilityDistribution{Univariate, Gaussian}) = ensureParameter!(dist, Val{:v}).params[:v] # unsafe variance
unsafeVar(dist::ProbabilityDistribution{Multivariate, Gaussian}) = diag(ensureParameter!(dist, Val{:v}).params[:v])

unsafeCov(dist::ProbabilityDistribution{Univariate, Gaussian}) = ensureParameter!(dist, Val{:v}).params[:v] # unsafe covariance
unsafeCov(dist::ProbabilityDistribution{Multivariate, Gaussian}) = deepcopy(ensureParameter!(dist, Val{:v}).params[:v])

function isProper(dist::ProbabilityDistribution{Univariate, Gaussian})
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

function isProper(dist::ProbabilityDistribution{Multivariate, Gaussian})
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

function isWellDefined(dist::ProbabilityDistribution{Univariate, Gaussian})
    # Check if dist is not underdetermined
    location_valid = isValid(dist, :m) || isValid(dist, :xi)
    scale_valid    = isValid(dist, :v) || isValid(dist, :w)

    return location_valid && scale_valid
end

function isWellDefined(dist::ProbabilityDistribution{Multivariate, Gaussian})
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

function sample(dist::ProbabilityDistribution{Univariate, Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :v))
    return sqrt(dist.params[:v])*randn() + dist.params[:m]
end

function sample(dist::ProbabilityDistribution{Multivariate, Gaussian})
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :v))
    return chol(dist.params[:v])' *randn(dims(dist)) + dist.params[:m]
end

function prod!( x::ProbabilityDistribution{Univariate, Gaussian},
                y::ProbabilityDistribution{Univariate, Gaussian},
                z::ProbabilityDistribution{Univariate, Gaussian}=ProbabilityDistribution(Univariate, Gaussian, xi=0.0, w=1.0))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    z.params[:m] = NaN
    z.params[:v] = NaN
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, Gaussian},
                            y::ProbabilityDistribution{Univariate, PointMass},
                            z::ProbabilityDistribution{Univariate, PointMass}=ProbabilityDistribution(Univariate, PointMass, m=0.0))

    z.params[:m] = y.params[:m]
    return z
end

function prod!( x::ProbabilityDistribution{Multivariate, Gaussian},
                y::ProbabilityDistribution{Multivariate, Gaussian},
                z::ProbabilityDistribution{Multivariate, Gaussian}=ProbabilityDistribution(Multivariate, Gaussian, xi=[NaN], w=[NaN].'))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    invalidate!(z, :m)
    invalidate!(z, :v)
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, Gaussian},
                            y::ProbabilityDistribution{Multivariate, PointMass},
                            z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN]))

    z.params[:m] = deepcopy(y.params[:m])
    return z
end

function ensureParameters!(dist::ProbabilityDistribution, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
function ensureParameter!{T<:VariateType}(dist::ProbabilityDistribution{T, Gaussian}, param::Type{Val{:m}})
    if !isValid(dist, :m)
        dist.params[:m] = ensureParameter!(dist, Val{:v}).params[:v] * dist.params[:xi]
    end
    return dist
end

function ensureParameter!{T<:VariateType}(dist::ProbabilityDistribution{T, Gaussian}, param::Type{Val{:v}})
    if !isValid(dist, :v)
        dist.params[:v] = cholinv(dist.params[:w])
    end
    return dist
end

function ensureParameter!{T<:VariateType}(dist::ProbabilityDistribution{T, Gaussian}, param::Type{Val{:xi}})
    if !isValid(dist, :xi)
        dist.params[:xi] = ensureParameter!(dist, Val{:w}).params[:w] * dist.params[:m]
    end
    return dist
end

function ensureParameter!{T<:VariateType}(dist::ProbabilityDistribution{T, Gaussian}, param::Type{Val{:w}})
    if !isValid(dist, :w)
        dist.params[:w] = cholinv(dist.params[:v])
    end
    return dist
end

function =={T<:VariateType}(t::ProbabilityDistribution{T, Gaussian}, u::ProbabilityDistribution{T, Gaussian})
    if t === u
        return true
    end
    ensureParameters!(t, (:xi, :w))
    ensureParameters!(u, (:xi, :w))
    return isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, Gaussian})
    return  0.5*log(unsafeCov(dist)) +
            0.5*log(2*pi) +
            0.5
end

function differentialEntropy(dist::ProbabilityDistribution{Multivariate, Gaussian})
    return  0.5*log(det(unsafeCov(dist))) +
            (dims(dist)/2)*log(2*pi) +
            (dims(dist)/2)
end