export SampleList

mutable struct SampleList <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    # For instant calculation of differential entropy
    diff_ent::Union{Number, Nothing}

    function SampleList(out, s, w; id=generateId(SampleList), diff_ent=nothing)
        @ensureVariables(out, s, w)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}(), diff_ent)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:s] = self.interfaces[2] = associate!(Interface(self), s)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{SampleList}) = "SampleList"

format(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = "$(slug(SampleList))(s=$(format(dist.params[:s])),w=$(format(dist.params[:w])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{SampleList}; s=[0.0], w=[1.0], diff_ent=nothing) = ProbabilityDistribution{Univariate, SampleList}(Dict{Symbol, Any}(:s=>s, :w=>w, :diff_ent=>diff_ent))
ProbabilityDistribution(::Type{Multivariate}, ::Type{SampleList}; s=[[0.0]], w=[1.0], diff_ent=nothing) = ProbabilityDistribution{Multivariate, SampleList}(Dict{Symbol, Any}(:s=>s, :w=>w, :diff_ent=>diff_ent))

dims(dist::ProbabilityDistribution{Univariate, SampleList}) = 1
dims(dist::ProbabilityDistribution{Multivariate, SampleList}) = length(dist.params[:s][1])

vague(::Type{SampleList}) = ProbabilityDistribution(Univariate, SampleList, s=rand(1000), w=ones(1000)/1000)

function vague(::Type{SampleList}, dims::Int64)
    s_list = Vector{Vector{Number}}(undef, 1000)
    for n=1:1000
        s_list[n] = rand(dims)
    end
    ProbabilityDistribution(Multivariate, SampleList, s=s_list, w=ones(1000)/1000)
end

unsafeMean(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = sum(log.(dist.params[:s]).*dist.params[:w])

# Unbiased (co)variance estimates
unsafeVar(dist::ProbabilityDistribution{Univariate, SampleList}) = (length(dist.params[:s])/(length(dist.params[:s])-1))*sum((dist.params[:s].-unsafeMean(dist)).^2 .*dist.params[:w])
unsafeCov(dist::ProbabilityDistribution{Univariate, SampleList}) = unsafeVar(dist)

unsafeVar(dist::ProbabilityDistribution{Multivariate, SampleList}) = diag(unsafeCov(dist))
function unsafeCov(dist::ProbabilityDistribution{Multivariate, SampleList})
    tot = zeros(dims(dist), dims(dist))
    m = unsafeMean(dist)
    samples = dist.params[:s]
    weights = dist.params[:w]
    for i=1:length(samples)
        tot += (samples[i].-m)*transpose(samples[i].-m) .* weights[i]
    end
    return (length(dist.params[:s])/(length(dist.params[:s])-1)).*tot
end

unsafeMeanCov(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = (unsafeMean(dist), unsafeCov(dist))

function unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, SampleList})
    all(0 .<= dist.params[:s] .< 1) || error("unsafeMirroredLogMean does not apply to variables outside of the range [0, 1)")

    return sum(log.(1 .- dist.params[:s]) .* dist.params[:w])
end

unsafeMeanVector(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

isProper(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = abs(sum(dist.params[:w]) - 1) < 0.001

# Prod functions are defined in such a way that bootstrap particle filter will be allowed
@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, SampleList},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    # Importance sampling - resampling
    log_pdf=(a) -> logPdf(x, a)
    w = exp.(log_pdf.(y.params[:s]))
    w = w .* y.params[:w]
    w = w./sum(w)
    # Compute effective number of particles to decide if resampling is needed
    n_eff = 1/sum(w.^2)
    if n_eff < length(w)/10
        weights = Weights(w)
        samples = sample(y.params[:s],weights,length(weights))

        z.params[:w] = ones(length(w))/length(w)
        z.params[:s] = samples
        return z
    else
        z.params[:w] = w
        z.params[:s] = y.params[:s]
        return z
    end
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate},
    y::ProbabilityDistribution{Multivariate, SampleList},
    z::ProbabilityDistribution{Multivariate, SampleList}=ProbabilityDistribution(Multivariate, SampleList, s=[[0.0]], w=[1.0]))

    # Importance sampling - resampling
    log_pdf=(a) -> logPdf(x, a)
    w = exp.(log_pdf.(y.params[:s]))
    w = w .* y.params[:w]
    w = w./sum(w)
    # Compute effective number of particles to decide if resampling is needed
    n_eff = 1/sum(w.^2)
    if n_eff < length(w)/10
        weights = Weights(w)
        samples = sample(y.params[:s],weights,length(weights))

        z.params[:w] = ones(length(w))/length(w)
        z.params[:s] = samples
        return z
    else
        z.params[:w] = w
        z.params[:s] = y.params[:s]
        return z
    end
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    sample_factor = Vector{Number}(undef, 1000)
    for i=1:1000
        sample_factor[i] = sample(x)
    end

    log_pdf=(a) -> y.params[:log_pdf](a)
    log_pdf_sf = log_pdf.(sample_factor)
    w = exp.(log_pdf_sf)
    H2 = log(sum(w)/1000) # To compute differential entropy
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    log_pdfx=(a) -> logPdf(x, a)
    H1 = -sum(w .* (log_pdfx.(sample_factor) .+ log_pdf_sf)) # To compute differential entropy
    z.params[:diff_ent] = H1+H2
    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Multivariate, SampleList}=ProbabilityDistribution(Multivariate, SampleList, s=[[0.0]], w=[1.0]))

    sample_factor = Vector{Vector{Number}}(undef, 1000)
    for i=1:1000
        sample_factor[i] = sample(x)
    end

    log_pdf=(a) -> y.params[:log_pdf](a)
    log_pdf_sf = log_pdf.(sample_factor)
    w = exp.(log_pdf_sf)
    H2 = log(sum(w)/1000) # To compute differential entropy
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    log_pdfx=(a) -> logPdf(x, a)
    H1 = -sum(w .* (log_pdfx.(sample_factor) .+ log_pdf_sf)) # To compute differential entropy
    z.params[:diff_ent] = H1+H2
    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Categorical},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))


    sample_factor = []
    for i=1:1000
        push!(sample_factor,sample(x))
    end
    
    log_pdf=(a) -> y.params[:log_pdf](a)
    log_pdf_sf = log_pdf.(sample_factor)
    w = exp.(log_pdf_sf)
    H2 = log(sum(w)/1000) # To compute differential entropy
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    log_pdfx=(a) -> logPdf(x, a)
    H1 = -sum(w .* (log_pdfx.(sample_factor) .+ log_pdf_sf)) # To compute differential entropy
    z.params[:diff_ent] = H1+H2
    return z
end

#This prod function is defined especially for nonconjugate inference
function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    sample_factor = Vector{Number}(undef, 1000)
    for i=1:1000
        sample_factor[i] = sample(x)
    end

    log_pdf=(a) -> logPdf(y,a)
    log_pdf_sf = log_pdf.(sample_factor)
    w = exp.(log_pdf_sf)
    H2 = log(sum(w)/1000) # To compute differential entropy
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    log_pdfx=(a) -> logPdf(x, a)
    H1 = -sum(w .* (log_pdfx.(sample_factor) .+ log_pdf_sf)) # To compute differential entropy
    z.params[:diff_ent] = H1+H2
    return z
end

# Differential entropy for SampleList
function differentialEntropy(dist::ProbabilityDistribution{V, SampleList} where V<:VariateType)
    dist.params[:diff_ent] === nothing && error("No applicable rule to approximate differential entropy for SampleList distribution")
    return dist.params[:diff_ent]
end
