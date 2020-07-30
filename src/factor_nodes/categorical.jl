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

format(dist::ProbabilityDistribution{Univariate, Categorical}) = "$(slug(Categorical))(p=$(format(dist.params[:p])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{Categorical}; p=[1/3, 1/3, 1/3]) = ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>p))
ProbabilityDistribution(::Type{Categorical}; p=[1/3, 1/3, 1/3]) = ProbabilityDistribution{Univariate, Categorical}(Dict(:p=>p))

dims(dist::ProbabilityDistribution{Univariate, Categorical}) = 1

vague(::Type{Categorical}, n_factors::Int64=3) = ProbabilityDistribution(Univariate, Categorical, p=(1/n_factors)*ones(n_factors))
vague(::Type{Categorical}, n_factors::Tuple) = vague(Categorical, n_factors[1])

isProper(dist::ProbabilityDistribution{Univariate, Categorical}) = (abs(sum(dist.params[:p])-1.) < 1e-6)

unsafeMean(dist::ProbabilityDistribution{Univariate, Categorical}) = deepcopy(dist.params[:p])

unsafeMeanVector(dist::ProbabilityDistribution{Univariate, Categorical}) = deepcopy(dist.params[:p])

function unsafeMode(dist::ProbabilityDistribution{Univariate, Categorical})
    i = findfirst(dist.params[:p] .== maximum(dist.params[:p])) # Index of first maximum
    m = zeros(length(dist.params[:p]))
    m[i] = 1.0

    return m
end

logPdf(dist::ProbabilityDistribution{Univariate, Categorical}, x) = sum(x .* log.(dist.params[:p]))

function sample(dist::ProbabilityDistribution{Univariate, Categorical})
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

function prod!( x::ProbabilityDistribution{Univariate, Categorical},
                y::ProbabilityDistribution{Univariate, Categorical},
                z::ProbabilityDistribution{Univariate, Categorical}=ProbabilityDistribution(Univariate, Categorical, p=ones(size(x.params[:p]))./length(x.params[:p])))

    # Multiplication of 2 categorical PMFs: p(z) = p(x) * p(y)
    z.params[:p][:] = clamp.(x.params[:p] .* y.params[:p], tiny, Inf) # Soften vanishing probabilities to enforce the convention [1, 0]*[0, 1] ∝ [1/2, 1/2]
    norm = sum(z.params[:p])
    z.params[:p] = z.params[:p]./norm

    return z
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Univariate, Categorical})
    -sum(dist.params[:p].*log.(clamp.(dist.params[:p], tiny, Inf))) # Soften vanishing propbabilities to enforce the convention 0 log 0 = 0
end

# Average energy functional
function averageEnergy(::Type{Categorical}, marg_out::ProbabilityDistribution, marg_p::ProbabilityDistribution{Multivariate})
    -sum(unsafeMean(marg_out).*unsafeLogMean(marg_p))
end
