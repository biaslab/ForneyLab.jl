export SampleList

mutable struct SampleList <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function SampleList(out, s, w; id=generateId(SampleList))
        @ensureVariables(out, s, w)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:s] = self.interfaces[2] = associate!(Interface(self), s)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{SampleList}) = "SampleList"

format(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = "$(slug(SampleList))(s=$(format(dist.params[:s])),w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{SampleList}; s=[0.0], w=[1.0]) = ProbabilityDistribution{Univariate, SampleList}(Dict{Symbol, Any}(:s=>s, :w=>w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{SampleList}; s=[[0.0]], w=[1.0]) = ProbabilityDistribution{Multivariate, SampleList}(Dict{Symbol, Any}(:s=>s, :w=>w))

dims(dist::ProbabilityDistribution{Univariate, SampleList}) = 1
dims(dist::ProbabilityDistribution{Multivariate, SampleList}) = length(dist.params[:s][1])

function vague(::Type{SampleList})
    n_samples = default_n_samples # Fixed number of samples

    return ProbabilityDistribution(Univariate, SampleList, s=rand(n_samples), w=ones(n_samples)/n_samples)
end

function vague(::Type{SampleList}, dims::Int64)
    n_samples = default_n_samples # Fixed number of samples

    s_list = Vector{Vector{Number}}(undef, n_samples)
    for n=1:n_samples
        s_list[n] = rand(dims)
    end
    
    return ProbabilityDistribution(Multivariate, SampleList, s=s_list, w=ones(n_samples)/n_samples)
end

unsafeMean(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = sum(log.(dist.params[:s]).*dist.params[:w])

# Unbiased (co)variance estimates
function unsafeVar(dist::ProbabilityDistribution{Univariate, SampleList})
    samples = dist.params[:s]
    n_samples = length(samples)

    return (n_samples/(n_samples - 1))*sum(dist.params[:w].*(samples .- unsafeMean(dist)).^2)
end

unsafeCov(dist::ProbabilityDistribution{Univariate, SampleList}) = unsafeVar(dist)

unsafeVar(dist::ProbabilityDistribution{Multivariate, SampleList}) = diag(unsafeCov(dist))
function unsafeCov(dist::ProbabilityDistribution{Multivariate, SampleList})
    samples = dist.params[:s]
    weights = dist.params[:w]

    n_samples = length(samples)
    m = unsafeMean(dist)
    tot = zeros(dims(dist), dims(dist))
    for i = 1:n_samples
        tot += (samples[i] .- m)*transpose(samples[i] .- m).*weights[i]
    end

    return (n_samples/(n_samples - 1)).*tot
end

unsafeMeanCov(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = (unsafeMean(dist), unsafeCov(dist))

function unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, SampleList})
    all(0 .<= dist.params[:s] .< 1) || error("unsafeMirroredLogMean does not apply to variables outside of the range [0, 1]")

    return sum(log.(1 .- dist.params[:s]) .* dist.params[:w])
end

unsafeMeanVector(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

isProper(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = abs(sum(dist.params[:w]) - 1) < 0.001

@symmetrical function prod!(
    x::ProbabilityDistribution{V}, # Includes function distributions
    y::ProbabilityDistribution{V, SampleList},
    z::ProbabilityDistribution{V, SampleList}=ProbabilityDistribution(V, SampleList, s=[0.0], w=[1.0])) where V<:VariateType

    samples = y.params[:s]
    n_samples = length(samples)
    log_samples_x = logPdf.([x], samples)

    # Compute sample weights
    w_raw_x = exp.(log_samples_x)
    w_prod = w_raw_x.*y.params[:w]
    weights = w_prod./sum(w_prod) # Normalize weights

    # Resample if required
    n_eff = 1/sum(weights.^2) # Effective number of particles
    if n_eff < n_samples/10
        samples = sample(samples, Weights(weights), n_samples) # Use StatsBase for resampling
        weights = ones(n_samples)./n_samples
    end

    # TODO: no entropy is computed here; include computation?
    z.params[:w] = weights
    z.params[:s] = samples

    return z
end

# Disambiguate beteen SampleList product and nonlinear product of any distribution with a Gaussian
# Following two definitions must be parameterized on separate Univariate and Multivariate types
@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, F},
    y::ProbabilityDistribution{Univariate, SampleList}) where F<:Gaussian

    z = ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0])

    return prod!(x, y, z) # Return a SampleList
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate, F},
    y::ProbabilityDistribution{Multivariate, SampleList}) where F<:Gaussian

    z = ProbabilityDistribution(Multivariate, SampleList, s=[[0.0]], w=[1.0])
    
    return prod!(x, y, z) # Return a SampleList
end

function sampleWeightsAndEntropy(x::ProbabilityDistribution, y::ProbabilityDistribution)
    n_samples = default_n_samples # Number of samples is fixed
    samples = sample(x, n_samples)

    # Apply log-pdf functions to the samples
    log_samples_x = logPdf.([x], samples)
    log_samples_y = logPdf.([y], samples)

    # Extract the sample weights
    w_raw = exp.(log_samples_y) # Unnormalized weights
    w_sum = sum(w_raw)
    weights = w_raw./w_sum # Normalize the raw weights

    # Compute the separate contributions to the entropy
    H_y = log(w_sum) - log(n_samples)
    H_x = -sum( weights.*(log_samples_x + log_samples_y) )
    entropy = H_x + H_y

    return (samples, weights, entropy)
end

# General product definition that returns a SampleList
function prod!(
    x::ProbabilityDistribution{V},
    y::ProbabilityDistribution{V},
    z::ProbabilityDistribution{V, SampleList} = ProbabilityDistribution(V, SampleList, s=[0.0], w=[1.0])) where V<:VariateType

    (samples, weights, entropy) = sampleWeightsAndEntropy(x, y)

    z.params[:s] = samples
    z.params[:w] = weights
    z.params[:entropy] = entropy

    return z
end

# Specialized product definition that accounts for the Categorical VariateType
@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Categorical},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Univariate, SampleList} = ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    (samples, weights, entropy) = sampleWeightsAndEntropy(x, y)

    z.params[:s] = samples
    z.params[:w] = weights
    z.params[:entropy] = entropy

    return z
end

# Differential entropy for SampleList
function differentialEntropy(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType
    haskey(dist.params, :entropy) || error("Missing entropy for SampleList; quantity is requested but not computed")
    
    return dist.params[:entropy] # Entropy is pre-computed during computation of the marginal
end
