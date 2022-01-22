export Dirichlet

"""
Description:

    Dirichlet factor node

    Multivariate:
    f(out, a) = Dir(out|a)
              = Γ(Σ_i a_i)/(Π_i Γ(a_i)) Π_i out_i^{a_i}
    where 'a' is a vector with every a_i > 0

    Matrix variate:
    f(out, a) = Π_k Dir(out|a_*k)
    where 'a' represents a left-stochastic matrix with every a_jk > 0

Interfaces:

    1. out
    2. a

Construction:

    Dirichlet(id=:some_id)
"""
mutable struct Dirichlet <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Dirichlet(out, a; id=generateId(Dirichlet))
        @ensureVariables(out, a)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:a] = self.interfaces[2] = associate!(Interface(self), a)

        return self
    end
end

slug(::Type{Dirichlet}) = "Dir"

format(dist::Distribution{V, Dirichlet}) where V<:VariateType = "$(slug(Dirichlet))(a=$(format(dist.params[:a])))"

Distribution(::Type{MatrixVariate}, ::Type{Dirichlet}; a=ones(3,3)) = Distribution{MatrixVariate, Dirichlet}(Dict(:a=>a))
Distribution(::Type{Multivariate}, ::Type{Dirichlet}; a=ones(3)) = Distribution{Multivariate, Dirichlet}(Dict(:a=>a))
Distribution(::Type{Dirichlet}; a=ones(3)) = Distribution{Multivariate, Dirichlet}(Dict(:a=>a))

dims(dist::Distribution{<:VariateType, Dirichlet}) = size(dist.params[:a])

vague(::Type{Dirichlet}, dims::Int64) = Distribution(Multivariate, Dirichlet, a=ones(dims))
vague(::Type{Dirichlet}, dims::Tuple{Int64}) = Distribution(Multivariate, Dirichlet, a=ones(dims))
vague(::Type{Dirichlet}, dims::Tuple{Int64, Int64}) = Distribution(MatrixVariate, Dirichlet, a=ones(dims))

isProper(dist::Distribution{V, Dirichlet}) where V<:VariateType = all(dist.params[:a] .> 0.0)

unsafeMean(dist::Distribution{Multivariate, Dirichlet}) = dist.params[:a]./sum(dist.params[:a])
unsafeMean(dist::Distribution{MatrixVariate, Dirichlet}) = dist.params[:a]./sum(dist.params[:a],dims=1) # Normalize columns

function unsafeLogMean(dist::Distribution{Multivariate, Dirichlet})
    a = clamp.(dist.params[:a], tiny, huge)
    return digamma.(a) .- digamma.(sum(a))
end
function unsafeLogMean(dist::Distribution{MatrixVariate, Dirichlet})
    a = clamp.(dist.params[:a], tiny, huge)
    return digamma.(a) .- digamma.(sum(a,dims=1)) # Normalize columns
end

logPdf(dist::Distribution{Multivariate, Dirichlet}, x) = sum((dist.params[:a].-1).*log.(x)) - sum(loggamma.(dist.params[:a])) + loggamma(sum(dist.params[:a]))
logPdf(dist::Distribution{MatrixVariate, Dirichlet}, x) = sum(sum((dist.params[:a].-1).*log.(x),dims=1) - sum(loggamma.(dist.params[:a]), dims=1) + loggamma.(sum(dist.params[:a],dims=1)))

function unsafeVar(dist::Distribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])
    return dist.params[:a].*(a_sum .- dist.params[:a])./(a_sum^2*(a_sum + 1.0))
end

function prod!( x::Distribution{V, Dirichlet},
                y::Distribution{V, Dirichlet},
                z::Distribution{V, Dirichlet}=Distribution(V, Dirichlet, a=ones(dims(x)))) where V<:VariateType

    z.params[:a] = x.params[:a] + y.params[:a] .- 1.0

    return z
end

@symmetrical function prod!(x::Distribution{Multivariate, Dirichlet},
                            y::Distribution{Multivariate, PointMass},
                            z::Distribution{Multivariate, PointMass}=Distribution(Multivariate, PointMass, m=[NaN]))

    a_y = clamp.(y.params[:m], 0.0, 1.0)
    z.params[:m] = deepcopy(a_y)

    return z
end

@symmetrical function prod!(x::Distribution{MatrixVariate, Dirichlet},
                            y::Distribution{MatrixVariate, PointMass},
                            z::Distribution{MatrixVariate, PointMass}=Distribution(MatrixVariate, PointMass, m=mat(NaN)))

    A_y = clamp.(y.params[:m], 0.0, 1.0)
    z.params[:m] = deepcopy(A_y)

    return z
end

function sample(dist::Distribution{Multivariate, Dirichlet})
    smpl = Vector{Float64}(undef, length(dist.params[:a]))
    for (i, alpha_i) in enumerate(dist.params[:a])
        smpl[i] = sample(Distribution(Univariate,Gamma,a=alpha_i,b=1.0))
    end
    smpl = smpl./sum(smpl)
    return smpl
end

naturalParams(dist::Distribution{Multivariate, Dirichlet}) = dist.params[:a] .- 1.0 # Variant 2 of https://en.wikipedia.org/wiki/Exponential_family
naturalParams(dist::Distribution{MatrixVariate, Dirichlet}) = vec(dist.params[:a]) .- 1.0

standardDistribution(::Type{Multivariate}, ::Type{Dirichlet}; η::Vector) = Distribution(Multivariate, Dirichlet, a=η.+1.0)
standardDistribution(::Type{MatrixVariate}, ::Type{Dirichlet}; η::Vector, dims::Tuple) = Distribution(MatrixVariate, Dirichlet, a=reshape(η, dims).+1.0) # Include dimensionality argument for rectangular case

logNormalizer(::Type{Multivariate}, ::Type{Dirichlet}; η::Vector) = sum(loggamma.(η.+1.0)) - loggamma(sum(η.+1.0))
logNormalizer(::Type{MatrixVariate}, ::Type{Dirichlet}; η::Vector, dims::Tuple) = sum(loggamma.(η.+1.0)) - sum(loggamma.(sum(reshape(η,dims).+1.0, dims=1)))

logPdf(V::Type{Multivariate}, F::Type{Dirichlet}, x::Vector; η::Vector) = log.(x)'*η - logNormalizer(V, F; η=η)
logPdf(V::Type{MatrixVariate}, F::Type{Dirichlet}, x::Matrix; η::Vector) = log.(vec(x))'*η - logNormalizer(V, F; η=η, dims=size(x))

# Entropy functional
function differentialEntropy(dist::Distribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])

    -sum( (dist.params[:a] .- 1.0).*(digamma.(dist.params[:a]) .- digamma.(a_sum)) ) -
    labsgamma(a_sum) +
    sum(labsgamma.(dist.params[:a]))
end

function differentialEntropy(dist::Distribution{MatrixVariate, Dirichlet})
    H = 0.0
    for k = 1:dims(dist)[2] # For all columns
        a_sum = sum(dist.params[:a][:,k])

        H += -sum( (dist.params[:a][:,k] .- 1.0).*(digamma.(dist.params[:a][:,k]) .- digamma.(a_sum)) ) -
        labsgamma(a_sum) +
        sum(labsgamma.(dist.params[:a][:,k]))
    end

    return H
end

# Average energy functional
function averageEnergy(::Type{Dirichlet}, marg_out::Distribution{Multivariate}, marg_a::Distribution{Multivariate, PointMass})
    a_sum = sum(marg_a.params[:m])

    -labsgamma(a_sum) +
    sum(labsgamma.(marg_a.params[:m])) -
    sum( (marg_a.params[:m] .- 1.0).*unsafeLogMean(marg_out) )
end

function averageEnergy(::Type{Dirichlet}, marg_out::Distribution{Multivariate}, marg_a::Distribution{Multivariate, SampleList})
    samples, weights = marg_a.params[:s], marg_a.params[:w]
    S = length(weights) #number of samples
    log_gamma_of_sum, sum_of_log_gamma = 0.0, 0.0
    
    for s=1:S
        log_gamma_of_sum += weights[s]*loggamma(sum(samples[s]))
        sum_of_log_gamma += weights[s]*sum(loggamma.(samples[s]))
    end

    return -log_gamma_of_sum + sum_of_log_gamma -
            sum((unsafeMeanVector(marg_a) .- 1.0).*unsafeLogMean(marg_out))
end

function averageEnergy(::Type{Dirichlet}, marg_out::Distribution{MatrixVariate}, marg_a::Distribution{MatrixVariate, PointMass})
    (dims(marg_out) == dims(marg_a)) || error("Distribution dimensions must agree")

    log_mean_marg_out = unsafeLogMean(marg_out)

    H = 0.0
    for k = 1:dims(marg_out)[2] # For all columns
        a_sum = sum(marg_a.params[:m][:,k])

        H += -labsgamma(a_sum) +
        sum(labsgamma.(marg_a.params[:m][:,k])) -
        sum( (marg_a.params[:m][:,k] .- 1.0).*log_mean_marg_out[:,k] )
    end

    return H
end
