export Abstract_dist

mutable struct Abstract_dist <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Abstract_dist(out, m, v, f; id=generateId(Abstract_dist))
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

slug(::Type{Abstract_dist}) = "Abstract_dist"

format(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = "$(slug(Abstract_dist))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])), f=$(format(dist.params[:f])))"

f_dummy(x) = x
ProbabilityDistribution(::Type{Univariate}, ::Type{Abstract_dist}; m=0.0, v=1.0, f=f_dummy) = ProbabilityDistribution{Univariate, Abstract_dist}(Dict(:m=>m, :v=>v, :f=>f))
ProbabilityDistribution(::Type{Munivariate}, ::Type{Abstract_dist}; m=[0.0, 0.0], v=[1.0 0;0 1.0], f=f_dummy) = ProbabilityDistribution{Univariate, Abstract_dist}(Dict(:m=>m, :v=>v, :f=>f))

dims(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = 1

samples(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = dist.params[:m] .+ sqrt(dist.params[:v]) .* randn(10000)
transformed_samples(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = dist.params[:f].(samples(dist))

unsafeMean(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = sum(transformed_samples(dist))/10000 # unsafe mean

unsafeLogMean(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = sum(log.(transformed_samples(dist)))/10000

function unsafeVar(dist::ProbabilityDistribution{Univariate, Abstract_dist})
    tr_samples = transformed_samples(dist)
    sum((tr_samples .- sum(tr_samples)/10000).^2)/9999
end

function unsafeMeanCov(dist::ProbabilityDistribution{Univariate, Abstract_dist})
    tr_samples = transformed_samples(dist)
    sum(tr_samples)/10000, sum((tr_samples .- sum(tr_samples)/10000).^2)/9999
end

function  unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, Abstract_dist})
    tr_samples = transformed_samples(dist)
    if all(0 .<= tr_samples .< 1)
        return sum(log.(1 .- transformed_samples(dist)))/10000
    else
        return error("Transformed samples are not in the range of [0,1]. Check your transformation function!")
    end
end

isProper(dist::ProbabilityDistribution{Univariate, Abstract_dist}) = (dist.params[:v] >= tiny)

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, Abstract_dist}=ProbabilityDistribution(Univariate, Abstract_dist, m=0.0, v=1.0, f=f_dummy)) where {F2<:Abstract_dist}

    y_m = y.params[:m]
    y_v = y.params[:v]
    y_f = y.params[:f]

    z.params[:m] = y_m
    z.params[:v] = y_v
    z.params[:f] = y_f

    return z
end
