abstract Gaussian <: SoftFactor

ProbabilityDistribution(::Type{Gaussian}) = ProbabilityDistribution(Gaussian, m=0.0, v=1.0)

function prod!( x::ProbabilityDistribution{Gaussian},
                y::ProbabilityDistribution{Gaussian},
                z::ProbabilityDistribution{Gaussian}=ProbabilityDistribution(Gaussian, xi=0.0, w=1.0))

    ensureParameters!(x, (:xi, :w))
    ensureParameters!(y, (:xi, :w))

    z.params[:m] = NaN
    z.params[:v] = NaN
    z.params[:xi] = x.params[:xi] + y.params[:xi]
    z.params[:w] = x.params[:w] + y.params[:w]

    return z
end

function ensureParameters!(dist::ProbabilityDistribution{Gaussian}, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:m}})
    if !haskey(dist.params, :m) || isnan(dist.params[:m])
        dist.params[:m] = ensureParameter!(dist, Val{:v}).params[:v] * dist.params[:xi]
    end
    return dist
end

function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:v}})
    if !haskey(dist.params, :v) || isnan(dist.params[:v])
        dist.params[:v] = 1/dist.params[:w]
    end
    return dist
end

function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:xi}})
    if !haskey(dist.params, :xi) || isnan(dist.params[:xi])
        dist.params[:xi] = ensureParameter!(dist, Val{:w}).params[:w] * dist.params[:m]
    end
    return dist
end

function ensureParameter!(dist::ProbabilityDistribution{Gaussian}, param::Type{Val{:w}})
    if !haskey(dist.params, :w) || isnan(dist.params[:w])
        dist.params[:w] = 1/dist.params[:v]
    end
    return dist
end

function ==(t::ProbabilityDistribution{Gaussian}, u::ProbabilityDistribution{Gaussian})
    ensureParameters!(t, (:xi, :w))
    ensureParameters!(u, (:xi, :w))
    return isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
end