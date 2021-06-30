export SetSampleList

"""
Description:

    SetProbDist is defined for messages that carry the posterior marginal
    information stored in a node to an edge and set its posterior marginal equal to
    the carried one. Required in SVI node.

"""
mutable struct SetSampleList <: SoftFactor
    params::Dict
end

slug(::Type{SetSampleList}) = "ssl"

format(dist::ProbabilityDistribution{V, SetSampleList}) where V<:VariateType = "$(dist.params)"

# Distribution constructors
ProbabilityDistribution(::Type{V}, ::Type{SetSampleList}; kwargs...) where V<:VariateType = ProbabilityDistribution{V, SetSampleList}(Dict{Symbol,Any}(kwargs))
ProbabilityDistribution(::Type{SetSampleList}; kwargs...) = ProbabilityDistribution{Univariate, SetSampleList}(Dict{Symbol,Any}(kwargs))


@symmetrical function prod!(x::ProbabilityDistribution{V, SetSampleList},
                            y::ProbabilityDistribution{V, F},
                            z::ProbabilityDistribution{V, SampleList} = ProbabilityDistribution(V, SampleList)) where {V<:VariateType, F<:FactorNode}

    thenode = currentGraph().nodes[x.params[:node_id]]
    thenode.message = y
    return deepcopy(x.params[:q])
end

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, SetSampleList},
                            y::ProbabilityDistribution{Univariate, F},
                            z::ProbabilityDistribution{Univariate, SampleList} = ProbabilityDistribution(Univariate, SampleList)) where {F<:Gaussian}

    thenode = currentGraph().nodes[x.params[:node_id]]
    thenode.message = y
    return deepcopy(x.params[:q])
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, SetSampleList},
                            y::ProbabilityDistribution{Multivariate, F},
                            z::ProbabilityDistribution{Multivariate, SampleList} = ProbabilityDistribution(Multivariate, SampleList)) where {F<:Gaussian}

    thenode = currentGraph().nodes[x.params[:node_id]]
    thenode.message = y
    return deepcopy(x.params[:q])
end
