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
ProbabilityDistribution(::Type{MatrixVariate}, ::Type{SampleList};s=[mat(0.0)], w=[1.0]) = ProbabilityDistribution{MatrixVariate, SampleList}(Dict{Symbol,Any}(:s=>s,:w=>w))

dims(dist::ProbabilityDistribution{Univariate, SampleList}) = 1
dims(dist::ProbabilityDistribution{Multivariate, SampleList}) = length(dist.params[:s][1])
dims(dist::ProbabilityDistribution{MatrixVariate, SampleList}) = size(dist.params[:s][1])

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

function vague(::Type{SampleList}, dims::Tuple{Int64,Int64})
    n_samples = default_n_samples # Fixed number of samples

    s_list = Vector{Matrix{Number}}(undef, n_samples)
    for n=1:n_samples
        s_list[n] = randn(dims)
    end

    return ProbabilityDistribution(Multivariate, SampleList, s=s_list, w=ones(n_samples)/n_samples)
end

unsafeMean(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = sum(log.(dist.params[:s]).*dist.params[:w])

unsafeMeanLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = sum(dist.params[:s].*log.(dist.params[:s]).*dist.params[:w])

function unsafeLogMean(dist::ProbabilityDistribution{Multivariate, SampleList})
    sum = zeros(length(dist.params[:s][1]))
    for i=1:length(dist.params[:s])
        sum = sum .+ log.(dist.params[:s][i]).*dist.params[:w][i]
    end
    return sum
end

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

function unsafeCov(dist::ProbabilityDistribution{MatrixVariate, SampleList})
    samples = dist.params[:s]
    weights = dist.params[:w]

    n_samples = length(samples)
    m = unsafeMean(dist)
    cov1 = zeros(dims(dist)[1],dims(dist)[1])
    cov2 = zeros(dims(dist)[2],dims(dist)[2])

    for i = 1:n_samples
        cov1 += ((samples[i] .- m))*transpose((samples[i] .- m)).*weights[i]
        cov2 += transpose((samples[i] .- m))*((samples[i] .- m)).*weights[i]
    end
    cov1 = (n_samples/(n_samples - 1)).*cov1
    cov2 = (n_samples/(n_samples - 1)).*cov2

    return kron(cov1, cov2)
end

unsafeMeanCov(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = (unsafeMean(dist), unsafeCov(dist))

function unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, SampleList})
    all(0 .<= dist.params[:s] .< 1) || error("unsafeMirroredLogMean does not apply to variables outside of the range [0, 1]")

    return sum(log.(1 .- dist.params[:s]) .* dist.params[:w])
end

unsafeMeanVector(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

isProper(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = abs(sum(dist.params[:w]) - 1) < 0.001

# prod of a pdf (or distribution) message and a SampleList message
# this function is capable to calculate entropy with SampleList messages in VMP setting
@symmetrical function prod!(
    x::ProbabilityDistribution{V}, # Includes function distributions
    y::ProbabilityDistribution{V, SampleList},
    z::ProbabilityDistribution{V, SampleList}=ProbabilityDistribution(V, SampleList, s=[0.0], w=[1.0])) where V<:VariateType

    # Suppose in the previous time step m1(pdf) and m2(pdf) messages collided.
    # The resulting collision m3 (sampleList) = m1*m2 is supposed to carry
    # the proposal (m1) and integrand (m2) distributions. m1 is the message from which
    # the samples are drawn. m2 is the message on which the samples are evaluated and
    # weights are calculated. In case Particle Filtering (BP), entropy will not be calculated
    # and in the first step there won't be any integrand information.
    if haskey(y.params, :logintegrand)
        # recall that we are calculating m3*m4. If m3 consists of integrand information
        # update it: new_integrand = m2*m3. This allows us to collide arbitrary number of beliefs
        # to approximate posterior and yet estimate the entropy.
        logIntegrand = (samples) -> y.params[:logintegrand](samples) .+ logPdf.([x], samples)
    else
        # If there is no integrand information before, set it to m4
        logIntegrand = (samples) -> logPdf.([x], samples)
    end

    samples = y.params[:s] # samples come from proposal (m1)
    n_samples = length(samples) # number of samples
    log_samples_x = logPdf.([x], samples) # evaluate samples in logm4, i.e. logm4(s)

    # Compute sample weights
    w_raw_x = exp.(log_samples_x) # m4(s)
    w_prod = w_raw_x.*y.params[:w] # update the weights of posterior w_unnormalized = m4(s)*w_prev
    weights = w_prod./sum(w_prod) # Normalize weights
    
    # Resample if required
    n_eff = 1/sum(weights.^2) # Effective number of particles
    if n_eff < n_samples/10
        samples = sample(samples, Weights(weights), n_samples) # Use StatsBase for resampling
        weights = ones(n_samples)./n_samples
    end

    # resulting posterior or message
    z.params[:w] = weights # set adjusted weights
    z.params[:s] = samples # samples are still coming from the same proposal
    z.params[:logintegrand] = logIntegrand # set integrand
    if haskey(y.params, :logproposal) && haskey(y.params, :unnormalizedweights)
        z.params[:unnormalizedweights] = w_raw_x.*y.params[:unnormalizedweights] # m4(s)*m2(s)
        logProposal = y.params[:logproposal] # m1
        z.params[:logproposal] = logProposal # m1
        # calculate entropy
        H_y = log(sum(w_raw_x.*y.params[:unnormalizedweights])) - log(n_samples) # log(sum_i(m4(s_i)*m2(s_i))/N)
        # -sum_i(w_i*log(m1(s_i)*m2(s_i)*m4(s_i)))
        H_x = -sum( weights.*(logProposal(samples) + log.(y.params[:unnormalizedweights]) + log_samples_x) )
        entropy = H_x + H_y
        z.params[:entropy] = entropy
    end

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

function sampleWeightsAndEntropy(x::ProbabilityDistribution{V,F}, y::ProbabilityDistribution) where {V<:VariateType, F<:Function}
    sampleWeightsAndEntropy(y, x)
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

    # Inform next step about the proposal and integrand to be used in entropy calculation in smoothing
    logproposal = (samples) -> logPdf.([x], samples)
    logintegrand = (samples) -> logPdf.([y], samples)

    return (samples, weights, w_raw, logproposal, logintegrand, entropy)
end

# General product definition that returns a SampleList
function prod!(
    x::ProbabilityDistribution{V},
    y::ProbabilityDistribution{V},
    z::ProbabilityDistribution{V, SampleList} = ProbabilityDistribution(V, SampleList, s=[0.0], w=[1.0])) where {V<:VariateType}

    (samples, weights, unnormalizedweights, logproposal, logintegrand, entropy) = sampleWeightsAndEntropy(x, y)

    z.params[:s] = samples
    z.params[:w] = weights
    z.params[:unnormalizedweights] = unnormalizedweights
    z.params[:logproposal] = logproposal
    z.params[:logintegrand] = logintegrand
    z.params[:entropy] = entropy

    return z
end

# Specialized product definition that accounts for the Categorical VariateType
@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Categorical},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Univariate, SampleList} = ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    (samples, weights, unnormalizedweights, logproposal, logintegrand, entropy) = sampleWeightsAndEntropy(x, y)

    z.params[:s] = samples
    z.params[:w] = weights
    z.params[:unnormalizedweights] = unnormalizedweights
    z.params[:logproposal] = logproposal
    z.params[:logintegrand] = logintegrand
    z.params[:entropy] = entropy

    return z
end

# Bootstrap samples
function bootstrap(dist_mean::ProbabilityDistribution{Univariate, SampleList}, dist_var::ProbabilityDistribution{Univariate, PointMass})
    s_m = dist_mean.params[:s] # Samples representing the mean
    N = length(s_m)
    v = dist_var.params[:m] # Fixed variance

    return sqrt(v)*randn(N) .+ s_m # New samples
end

function bootstrap(dist_mean::ProbabilityDistribution{Multivariate, SampleList}, dist_var::ProbabilityDistribution{MatrixVariate, PointMass})
    d = dims(dist_mean)
    s_m = dist_mean.params[:s] # Samples representing the mean
    N = length(s_m)
    V = dist_var.params[:m] # Fixed variance
    U = (cholesky(default_cholesky_mode, V)).U # Precompute Cholesky

    return [U' *randn(d) + s_m[i] for i in 1:N] # New samples
end

function bootstrap(dist_mean::ProbabilityDistribution{Univariate, <:Gaussian}, dist_var::ProbabilityDistribution{Univariate, SampleList})
    s_v = dist_var.params[:s] # Samples representing the variance
    N = length(s_v)
    (m, v) = unsafeMeanCov(dist_mean)
    s_u = sqrt.(s_v .+ v) # Standard deviation for each variance sample

    return s_u.*randn(N) .+ m # New samples
end

function bootstrap(dist_mean::ProbabilityDistribution{Multivariate, <:Gaussian}, dist_var::ProbabilityDistribution{MatrixVariate, SampleList})
    d = dims(dist_mean)
    s_V = dist_var.params[:s] # Samples representing the covariance
    N = length(s_V)
    (m, V) = unsafeMeanCov(dist_mean)
    s_U = [(cholesky(default_cholesky_mode, s_V[i] + V)).U for i in 1:N] # Precompute Cholesky for each covariance sample; this can be expensive

    return [s_U[i]' *randn(d) + m for i in 1:N] # New samples
end

#Sampling from a distribution in ForneyLab returns equally weighted samples from the distribution
#To retain the unified standard procedure, we allow sampling from SampleList not through directly returning
#sample and weight parameters but drawing samples from the existing list of samples according to weights.
# Inverse-transform sampling
sample(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sample(dist.params[:s], Weights(dist.params[:w]))

# Differential entropy for SampleList
function differentialEntropy(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType
    haskey(dist.params, :entropy) || error("Missing entropy for SampleList; quantity is requested but not computed")

    return dist.params[:entropy] # Entropy is pre-computed during computation of the marginal
end
