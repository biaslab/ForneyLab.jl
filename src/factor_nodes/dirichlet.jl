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

format(dist::ProbabilityDistribution{V, Dirichlet}) where V<:VariateType = "$(slug(Dirichlet))(a=$(format(dist.params[:a])))"

ProbabilityDistribution(::Type{MatrixVariate}, ::Type{Dirichlet}; a=ones(3,3)) = ProbabilityDistribution{MatrixVariate, Dirichlet}(Dict(:a=>a))
ProbabilityDistribution(::Type{Multivariate}, ::Type{Dirichlet}; a=ones(3)) = ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>a))
ProbabilityDistribution(::Type{Dirichlet}; a=ones(3)) = ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>a))

dims(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = length(dist.params[:a])
dims(dist::ProbabilityDistribution{MatrixVariate, Dirichlet}) = size(dist.params[:a])

vague(::Type{Dirichlet}, dims::Int64) = ProbabilityDistribution(Multivariate, Dirichlet, a=ones(dims))
vague(::Type{Dirichlet}, dims::Tuple{Int64}) = ProbabilityDistribution(Multivariate, Dirichlet, a=ones(dims))
vague(::Type{Dirichlet}, dims::Tuple{Int64, Int64}) = ProbabilityDistribution(MatrixVariate, Dirichlet, a=ones(dims))

isProper(dist::ProbabilityDistribution{V, Dirichlet}) where V<:VariateType = all(dist.params[:a] .> 0.0)

unsafeMean(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = dist.params[:a]./sum(dist.params[:a])
unsafeMean(dist::ProbabilityDistribution{MatrixVariate, Dirichlet}) = dist.params[:a]./sum(dist.params[:a],dims=1) # Normalize columns

unsafeLogMean(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = digamma.(dist.params[:a]) .- digamma.(sum(dist.params[:a]))
unsafeLogMean(dist::ProbabilityDistribution{MatrixVariate, Dirichlet}) = digamma.(dist.params[:a]) .- digamma.(sum(dist.params[:a],dims=1)) # Normalize columns

logPdf(dist::ProbabilityDistribution{Multivariate, Dirichlet}, x) = sum((dist.params[:a].-1).*log.(x)) - sum(loggamma.(dist.params[:a])) + loggamma(sum(dist.params[:a]))
logPdf(dist::ProbabilityDistribution{MatrixVariate, Dirichlet}, x) = sum(sum((dist.params[:a].-1).*log.(x),dims=1) - sum(loggamma.(dist.params[:a]), dims=1) + loggamma.(sum(dist.params[:a],dims=1)))

function unsafeVar(dist::ProbabilityDistribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])
    return dist.params[:a].*(a_sum .- dist.params[:a])./(a_sum^2*(a_sum + 1.0))
end

function prod!( x::ProbabilityDistribution{V, Dirichlet},
                y::ProbabilityDistribution{V, Dirichlet},
                z::ProbabilityDistribution{V, Dirichlet}=ProbabilityDistribution(V, Dirichlet, a=ones(dims(x)))) where V<:VariateType

    z.params[:a] = x.params[:a] + y.params[:a] .- 1.0

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, Dirichlet},
                            y::ProbabilityDistribution{Multivariate, PointMass},
                            z::ProbabilityDistribution{Multivariate, PointMass}=ProbabilityDistribution(Multivariate, PointMass, m=[NaN]))

    all(0.0 .<= y.params[:m] .<= 1.0) || error("PointMass location entries $(y.params[:m]) should all be between 0 and 1")
    isapprox(sum(y.params[:m]), 1.0) || error("Pointmass location entries $(y.params[:m]) should sum to one")
    z.params[:m] = deepcopy(y.params[:m])

    return z
end

@symmetrical function prod!(x::ProbabilityDistribution{MatrixVariate, Dirichlet},
                            y::ProbabilityDistribution{MatrixVariate, PointMass},
                            z::ProbabilityDistribution{MatrixVariate, PointMass}=ProbabilityDistribution(MatrixVariate, PointMass, m=mat(NaN)))

    all(0.0 .<= y.params[:m] .<= 1.0) || error("PointMass location entries $(y.params[:m]) should all be between 0 and 1")
    for k = 1:dims(y)[2] # For all columns
        isapprox(sum(y.params[:m][:,k]), 1.0) || error("Pointmass location entries $(y.params[:m][:,k]) of column $k should sum to one")
    end

    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function sample(dist::ProbabilityDistribution{Multivariate, Dirichlet})
    smpl = Vector{Float64}(undef, length(dist.params[:a]))
    for (i, alpha_i) in enumerate(dist.params[:a])
        smpl[i] = sample(ProbabilityDistribution(Univariate,Gamma,a=alpha_i,b=1.0))
    end
    smpl = smpl./sum(smpl)
    return smpl
end

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])

    -sum( (dist.params[:a] .- 1.0).*(digamma.(dist.params[:a]) .- digamma.(a_sum)) ) -
    labsgamma(a_sum) +
    sum(labsgamma.(dist.params[:a]))
end

function differentialEntropy(dist::ProbabilityDistribution{MatrixVariate, Dirichlet})
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
function averageEnergy(::Type{Dirichlet}, marg_out::ProbabilityDistribution{Multivariate}, marg_a::ProbabilityDistribution{Multivariate, PointMass})
    a_sum = sum(marg_a.params[:m])

    -labsgamma(a_sum) +
    sum(labsgamma.(marg_a.params[:m])) -
    sum( (marg_a.params[:m] .- 1.0).*unsafeLogMean(marg_out) )
end

function averageEnergy(::Type{Dirichlet}, marg_out::ProbabilityDistribution{Multivariate}, marg_a::ProbabilityDistribution{Multivariate, SampleList})
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

function averageEnergy(::Type{Dirichlet}, marg_out::ProbabilityDistribution{MatrixVariate}, marg_a::ProbabilityDistribution{MatrixVariate, PointMass})
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
