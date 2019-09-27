export RGMP_dist

mutable struct RGMP_dist <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function RGMP_dist(out, m, v, f; id=generateId(RGMP_dist))
        @ensureVariables(out, m, v, f)
        self = new(id, Array{Interface}(undef, 4), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
        self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)
        self.i[:f] = self.interfaces[4] = associate!(Interface(self), f)

        return self
    end
end

slug(::Type{RGMP_dist}) = "RGMP_dist"

format(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = "$(slug(RGMP_dist))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])), f=$(format(dist.params[:f])))"

f_dummy(x) = x
ProbabilityDistribution(::Type{Univariate}, ::Type{RGMP_dist}; m=0.0, v=1.0, f=f_dummy) = ProbabilityDistribution{Univariate, RGMP_dist}(Dict(:m=>m, :v=>v, :f=>f))
ProbabilityDistribution(::Type{Munivariate}, ::Type{RGMP_dist}; m=[0.0, 0.0], v=[1.0 0;0 1.0], f=f_dummy) = ProbabilityDistribution{Univariate, RGMP_dist}(Dict(:m=>m, :v=>v, :f=>f))

dims(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = 1

samples(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = dist.params[:m] .+ sqrt(dist.params[:v]) .* randn(100)
transformed_samples(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = dist.params[:f].(samples(dist))

unsafeMean(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = sum(transformed_samples(dist))/100 # unsafe mean

unsafeLogMean(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = sum(log.(transformed_samples(dist)))/100

function unsafeVar(dist::ProbabilityDistribution{Univariate, RGMP_dist})
    tr_samples = transformed_samples(dist)
    sum((tr_samples .- sum(tr_samples)/100).^2)/99
end

function unsafeMeanCov(dist::ProbabilityDistribution{Univariate, RGMP_dist})
    tr_samples = transformed_samples(dist)
    sum(tr_samples)/100, sum((tr_samples .- sum(tr_samples)/100).^2)/99
end

function  unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, RGMP_dist})
    tr_samples = transformed_samples(dist)
    if tr_samples[0 .< tr_samples .< 1]
        return sum(log.(1 .- transformed_samples(dist)))/100
    else
        return error("Transformed samples are not in the range of [0,1]. Check your transformation function!")
    end
end

isProper(dist::ProbabilityDistribution{Univariate, RGMP_dist}) = (dist.params[:v] >= tiny)

function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, RGMP_dist}=ProbabilityDistribution(Univariate, RGMP_dist, m=0.0, v=1.0, f=f_dummy)) where {F2<:RGMP_dist}

    y_m = y.params[:m]
    y_v = y.params[:v]
    y_f = y.params[:f]

    z.params[:m] = y_m
    z.params[:v] = y_v
    z.params[:f] = y_f

    return z
end

function prod!(
    y::ProbabilityDistribution{Univariate, F1},
    x::ProbabilityDistribution{Univariate},
    z::ProbabilityDistribution{Univariate, RGMP_dist}=ProbabilityDistribution(Univariate, RGMP_dist, m=0.0, v=1.0, f=f_dummy)) where {F1<:RGMP_dist}

    y_m = y.params[:m]
    y_v = y.params[:v]
    y_f = y.params[:f]

    z.params[:m] = y_m
    z.params[:v] = y_v
    z.params[:f] = y_f

    return z
end
