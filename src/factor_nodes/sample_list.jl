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

ProbabilityDistribution(::Type{Univariate}, ::Type{SampleList}; s=[0.0], w=[1.0]) = ProbabilityDistribution{Univariate, SampleList}(Dict(:s=>s, :w=>w))
ProbabilityDistribution(::Type{Multivariate}, ::Type{SampleList}; s=[[0.0]], w=[1.0]) = ProbabilityDistribution{Multivariate, SampleList}(Dict(:s=>s, :w=>w))

dims(dist::ProbabilityDistribution{Univariate, SampleList}) = 1
dims(dist::ProbabilityDistribution{Multivariate, SampleList}) = length(dist.params[:s][1])

vague(::Type{SampleList}) = ProbabilityDistribution(Univariate, SampleList, s=rand(1000), w=ones(1000)/1000)
function vague(::Type{SampleList}, dims::Int64)
    s_list = []
    for n=1:1000
        append!(s_list,rand(dims))
    end
    ProbabilityDistribution(Multivariate, SampleList, s=s_list, w=ones(1000)/1000)
end

unsafeMean(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = sum(dist.params[:s].*dist.params[:w])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = sum(log.(dist.params[:s]).*dist.params[:w])

#Unbiased variance estimates
unsafeVar(dist::ProbabilityDistribution{Univariate, SampleList}) = (length(dist.params[:s])/(length(dist.params[:s])-1))*sum((dist.params[:s].-unsafeMean(dist)).^2 .*dist.params[:w])

function unsafeVar(dist::ProbabilityDistribution{Multivariate, SampleList})
    tot = zeros(dims(dist), dims(dist))
    m = unsafeMean(dist)
    samples = dist.params[:s]
    weights = dist.params[:w]
    for i=1:length(samples)
        tot += (samples[i].-m)*transpose(samples[i].-m) .* weights[i]
    end
    return (length(dist.params[:s])/(length(dist.params[:s])-1)).*tot
end

unsafeMeanCov(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = (unsafeMean(dist), unsafeVar(dist))

function unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, SampleList})
    all(0 .<= dist.params[:s] .< 1) || error("unsafeMirroredLogMean does not apply to variables outside of the range [0, 1)")

    return sum(log.(1 .- dist.params[:s]) .* dist.params[:w])
end

isProper(dist::ProbabilityDistribution{V, SampleList}) where V<:VariateType = abs(sum(dist.params[:w]) - 1) < 0.001

#Prod functions are defined in such a way that bootstrap particle filter will be allowed
@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, SampleList},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0], w=[1.0]))

    #Importance sampling - resampling
    log_pdf=(a) -> logPdf(x, a)
    w = exp.(log_pdf.(y.params[:s]))
    w = w .* y.params[:w]
    w = w./sum(w)
    #compute effective number of particles to decide if resampling is needed
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

    #Importance sampling - resampling
    log_pdf=(a) -> logPdf(x, a)
    w = exp.(log_pdf.(y.params[:s]))
    w = w .* y.params[:w]
    w = w./sum(w)
    #compute effective number of particles to decide if resampling is needed
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

    sample_factor = []
    for i=1:1000
        push!(sample_factor,sample(x))
    end

    log_pdf=(a) -> y.params[:log_pdf](a)
    w = exp.(log_pdf.(sample_factor))
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    return z
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Multivariate},
    y::ProbabilityDistribution{Multivariate, Function},
    z::ProbabilityDistribution{Multivariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[[0.0]], w=[1.0]))

    sample_factor = []
    for i=1:1000
        push!(sample_factor,sample(x))
    end

    log_pdf=(a) -> y.params[:log_pdf](a)
    w = exp.(log_pdf.(sample_factor))
    w = w./sum(w)
    z.params[:w] = w
    z.params[:s] = sample_factor
    return z
end
