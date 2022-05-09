export Categorical

"""
Description:

    Categorical factor node

    The categorical node defines a one-dimensional probability
    distribution over the normal basis vectors of dimension d

    out ∈ {0, 1}^d where Σ_k out_k = 1
    p ∈ [0, 1]^d, where Σ_k p_k = 1

    f(out, p) = Cat(out | p)
              = Π_i p_i^{out_i}

Interfaces:

    1. out
    2. p

Construction:

    Categorical(id=:some_id)
"""
mutable struct Categorical <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Categorical(out, p; id=generateId(Categorical))
        @ensureVariables(out, p)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:p] = self.interfaces[2] = associate!(Interface(self), p)

        return self
    end
end

slug(::Type{Categorical}) = "Cat"

format(dist::Distribution{Univariate, Categorical}) = "$(slug(Categorical))(p=$(format(dist.params[:p])))"

Distribution(::Type{Univariate}, ::Type{Categorical}; p=[1/3, 1/3, 1/3]) = Distribution{Univariate, Categorical}(Dict(:p=>p))
Distribution(::Type{Categorical}; p=[1/3, 1/3, 1/3]) = Distribution{Univariate, Categorical}(Dict(:p=>p))

dims(dist::Distribution{Univariate, Categorical}) = ()

vague(::Type{Categorical}, n_factors::Int64=3) = Distribution(Univariate, Categorical, p=(1/n_factors)*ones(n_factors))
vague(::Type{Categorical}, n_factors::Tuple) = vague(Categorical, n_factors[1])

isProper(dist::Distribution{Univariate, Categorical}) = (abs(sum(dist.params[:p])-1.) < 1e-6)

unsafeMean(dist::Distribution{Univariate, Categorical}) = deepcopy(dist.params[:p])

unsafeMeanVector(dist::Distribution{Univariate, Categorical}) = deepcopy(dist.params[:p])

function unsafeMode(dist::Distribution{Univariate, Categorical})
    i = findfirst(dist.params[:p] .== maximum(dist.params[:p])) # Index of first maximum
    m = zeros(length(dist.params[:p]))
    m[i] = 1.0

    return m
end

logPdf(dist::Distribution{Univariate, Categorical}, x) = sum(x .* log.(dist.params[:p]))

function sample(dist::Distribution{Univariate, Categorical})
    isProper(dist) || error("Cannot sample from improper distribution $dist")

    p = dist.params[:p]
    p_sum = cumsum(p);
    z = rand()

    d = length(p)
    idx = 0
    for i = 1:d
        if z < p_sum[i]
            idx = i
            break
        end
        idx = d
    end

    x = spzeros(Float64, d)
    x[idx] = 1.0

    return x
end

naturalParams(dist::Distribution{Univariate, Categorical}) = vcat(log.(dist.params[:p][1:end-1]./dist.params[:p][end]), 0.0) # Variant 3 of https://en.wikipedia.org/wiki/Exponential_family

standardDistribution(V::Type{Univariate}, F::Type{Categorical}; η::Vector) = Distribution(V, F, p=exp.(η)./sum(exp.(η)))

logNormalizer(::Type{Univariate}, ::Type{Categorical}; η::Vector) = log(sum(exp.(η)))

logPdf(V::Type{Univariate}, F::Type{Categorical}, x::AbstractVector; η::Vector) = x'*η - logNormalizer(V, F; η=η)

function prod!( x::Distribution{Univariate, Categorical},
                y::Distribution{Univariate, Categorical},
                z::Distribution{Univariate, Categorical}=Distribution(Univariate, Categorical, p=ones(size(x.params[:p]))./length(x.params[:p])))

    # Multiplication of 2 categorical PMFs: p(z) = p(x) * p(y)
    z.params[:p][:] = clamp.(x.params[:p] .* y.params[:p], tiny, Inf) # Soften vanishing probabilities to enforce the convention [1, 0]*[0, 1] ∝ [1/2, 1/2]
    norm = sum(z.params[:p])
    z.params[:p] = z.params[:p]./norm

    return z
end

@symmetrical function prod!(x::Distribution{Univariate, Categorical},
                            y::Distribution{Multivariate, PointMass},
                            z::Distribution{Multivariate, PointMass}=Distribution(Multivariate, PointMass, m=[0.0]))

    z.params[:m] = deepcopy(y.params[:m])

    return z
end


# Entropy functional
function differentialEntropy(dist::Distribution{Univariate, Categorical})
    -sum(dist.params[:p].*log.(clamp.(dist.params[:p], tiny, Inf))) # Soften vanishing propbabilities to enforce the convention 0 log 0 = 0
end

# Average energy functional
function averageEnergy(::Type{Categorical}, marg_out::Distribution, marg_p::Distribution{Multivariate})
    -sum(unsafeMean(marg_out).*unsafeLogMean(marg_p))
end

function averageEnergy(::Type{Categorical}, marg_out::Distribution{Multivariate,PointMass}, marg_p::Distribution{Multivariate})
    -sum(marg_out.params[:m].*unsafeLogMean(marg_p))
end
