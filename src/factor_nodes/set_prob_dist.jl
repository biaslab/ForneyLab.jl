export SetProbDist

"""
Description:

    SetProbDist is defined for messages that carry the posterior marginal
    informations stored in a node to an edge. Required in SVI node.

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
    thenode = currentGraph().nodes[x.params[:node_id]]
    marg = deepcopy(thenode.q_memory)
    thenode.q_memory = thenode.q
    return marg
end

# For some reason Gaussian messages do not enter in F<:FactorNode
@symmetrical function prod!(x::ProbabilityDistribution{Univariate, SetProbDist},
                            y::ProbabilityDistribution{Univariate, F},
                            z::ProbabilityDistribution{Univariate, F} = ProbabilityDistribution(Univariate,F)) where {F<:Gaussian}
    thenode = currentGraph().nodes[x.params[:node_id]]
    marg = deepcopy(thenode.q_memory)
    thenode.q_memory = thenode.q
    return marg
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, SetProbDist},
                            y::ProbabilityDistribution{Multivariate, F},
                            z::ProbabilityDistribution{Multivariate, F} = ProbabilityDistribution(Multivariate,F)) where {F<:Gaussian}
    thenode = currentGraph().nodes[x.params[:node_id]]
    marg = deepcopy(thenode.q_memory)
    thenode.q_memory = thenode.q
    return marg
end
