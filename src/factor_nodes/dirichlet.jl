export Dirichlet

"""
Description:
    Dirichlet factor node

    Real vector
    a .> 0
    
    f(out, a) = Dirichlet(out|a)

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

    function Dirichlet(out::Variable, a::Variable; id=generateId(Dirichlet))
        self = new(id, Array{Interface}(2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:a] = self.interfaces[2] = associate!(Interface(self), a)

        return self
    end
end

slug(::Type{Dirichlet}) = "Dir"

format(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = "$(slug(Dirichlet))(a=$(format(dist.params[:a])))"

ProbabilityDistribution(::Type{Multivariate}, ::Type{Dirichlet}; a=ones(3)) = ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>a))
ProbabilityDistribution(::Type{Dirichlet}; a=ones(3)) = ProbabilityDistribution{Multivariate, Dirichlet}(Dict(:a=>a))

dims(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = length(dist.params[:a])

vague(::Type{Dirichlet}, dims::Int64) = ProbabilityDistribution(Multivariate, Dirichlet, a=ones(dims))

isProper(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = all(dist.params[:a] .> 0.0)

unsafeMean(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = dist.params[:a]./sum(dist.params[:a])

unsafeLogMean(dist::ProbabilityDistribution{Multivariate, Dirichlet}) = digamma(dist.params[:a]) - digamma(sum(dist.params[:a]))

function unsafeVar(dist::ProbabilityDistribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])
    return dist.params[:a].*(a_sum - dist.params[:a])./(a_sum^2*(a_sum + 1.0))
end

function prod!( x::ProbabilityDistribution{Multivariate, Dirichlet},
                y::ProbabilityDistribution{Multivariate, Dirichlet},
                z::ProbabilityDistribution{Multivariate, Dirichlet}=ProbabilityDistribution(Multivariate, Dirichlet, a=ones(dims(x))))

    z.params[:a] = x.params[:a] + y.params[:a] - 1.0

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

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Multivariate, Dirichlet})
    a_sum = sum(dist.params[:a])

    -sum( (dist.params[:a] - 1.0).*(digamma(dist.params[:a]) - digamma(a_sum)) ) -
    lgamma(a_sum) +
    sum( lgamma(dist.params[:a]) )
end

# Average energy functional
function averageEnergy(::Type{Dirichlet}, marg_out::ProbabilityDistribution{Multivariate}, marg_a::ProbabilityDistribution{Multivariate, PointMass})
    a_sum = sum(marg_a.params[:m])

    -lgamma(a_sum) +
    sum( lgamma(marg_a.params[:m]) ) -
    sum( (marg_a.params[:m] - 1.0).*unsafeLogMean(marg_out) )
end
