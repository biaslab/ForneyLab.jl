export SampleList

mutable struct SampleList <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function SampleList(out, s; id=generateId(SampleList))
        @ensureVariables(out, s)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:s] = self.interfaces[2] = associate!(Interface(self), s)

        return self
    end
end

slug(::Type{SampleList}) = "SampleList"

format(dist::ProbabilityDistribution{Univariate, SampleList}) = "$(slug(SampleList))(s=$(format(dist.params[:s])))"

ProbabilityDistribution(::Type{Univariate}, ::Type{SampleList}; s=[0.0]) = ProbabilityDistribution{Univariate, SampleList}(Dict(:s=>s))

dims(dist::ProbabilityDistribution{Univariate, SampleList}) = 1

unsafeMean(dist::ProbabilityDistribution{Univariate, SampleList}) = mean(dist.params[:s])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, SampleList}) = mean(log.(dist.params[:s]))

unsafeVar(dist::ProbabilityDistribution{Univariate, SampleList}) = var(dist.params[:s]) # Compute corrected variance

unsafeMeanCov(dist::ProbabilityDistribution{Univariate, SampleList}) = (mean(dist.params[:s]), var(dist.params[:s]))

function unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, SampleList})
    all(0 .<= dist.params[:s] .< 1) || error("unsafeMirroredLogMean does not apply to variables outside of the range [0, 1)")

    return mean(log.(1 .- dist.params[:s]))
end

isProper(dist::ProbabilityDistribution{Univariate, SampleList}) = true

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, SampleList},
    z::ProbabilityDistribution{Univariate, SampleList}=ProbabilityDistribution(Univariate, SampleList, s=[0.0]))

    #Importance sampling - resampling
    log_pdf=(a) -> logPdf(x, a)
    weights = exp.(log_pdf.(y.params[:s]))
    weights = Weights(weights./sum(weights))
    samples = sample(y.params[:s],weights,length(weights))

    z.params[:s] = samples
    return z
end
