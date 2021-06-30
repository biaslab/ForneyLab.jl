export SetProbDist

"""
Description:

    SetProbDist is defined for messages that carry the posterior marginal
    information stored in a node to an edge and set its posterior marginal equal to
    the carried one. Required in SVI node.

"""
mutable struct SetProbDist <: SoftFactor
    params::Dict
end

slug(::Type{SetProbDist}) = "spd"

format(dist::ProbabilityDistribution{V, SetProbDist}) where V<:VariateType = "$(dist.params)"

# Distribution constructors
ProbabilityDistribution(::Type{V}, ::Type{SetProbDist}; kwargs...) where V<:VariateType = ProbabilityDistribution{V, SetProbDist}(Dict{Symbol,Any}(kwargs))
ProbabilityDistribution(::Type{SetProbDist}; kwargs...) = ProbabilityDistribution{Univariate, SetProbDist}(Dict{Symbol,Any}(kwargs))


@symmetrical function prod!(x::ProbabilityDistribution{V, SetProbDist},
                            y::ProbabilityDistribution{V, F},
                            z::ProbabilityDistribution{V, F} = ProbabilityDistribution(V,F)) where {V<:VariateType, F<:FactorNode}
    return deepcopy(x.params[:q])
end

# For some reason Gaussian messages do not enter in F<:FactorNode
@symmetrical function prod!(x::ProbabilityDistribution{Univariate, SetProbDist},
                            y::ProbabilityDistribution{Univariate, F},
                            z::ProbabilityDistribution{Univariate, F} = ProbabilityDistribution(Univariate,F)) where {F<:Gaussian}
    return deepcopy(x.params[:q])
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, SetProbDist},
                            y::ProbabilityDistribution{Multivariate, F},
                            z::ProbabilityDistribution{Multivariate, F} = ProbabilityDistribution(Univariate,F)) where {F<:Gaussian}
    return deepcopy(x.params[:q])
end




# @symmetrical function prod!(x::ProbabilityDistribution{Univariate, SetProbDist},
#                             y::ProbabilityDistribution{Univariate, F1},
#                             z::ProbabilityDistribution{Univariate, SampleList} = ProbabilityDistribution(Univariate, SampleList)) where {F1<:Gaussian}
#     @show x.params[:q]
#     return deepcopy(x.params[:q])
# end
