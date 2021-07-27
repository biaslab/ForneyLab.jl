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


# isUnivariateGaussian(dist::ProbabilityDistribution{Univariate, F}) where F<:Gaussian = true
# isUnivariateGaussian(dist::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:FactorFunction} = false
# isMultivariateGaussian(dist::ProbabilityDistribution{Multivariate, F}) where F<:Gaussian = true
# isMultivariateGaussian(dist::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:FactorFunction} = false

@symmetrical function prod!(x::ProbabilityDistribution{V, SetSampleList},
                            y::ProbabilityDistribution{V, F},
                            z::ProbabilityDistribution{V, SampleList} = ProbabilityDistribution(V, SampleList)) where {V<:VariateType, F<:FactorFunction}

    thenode = currentGraph().nodes[x.params[:node_id]]

    samples_ = [sample(prob, thenode.num_samples) for prob in thenode.q_memory]
    samples = thenode.g.(samples_...)
    weights = ones(thenode.num_samples)/thenode.num_samples
    thenode.q_memory = thenode.q
    return ProbabilityDistribution(V, SampleList, s=samples, w=weights)
end

# @symmetrical function prod!(x::ProbabilityDistribution{V, SetSampleList},
#                             y::ProbabilityDistribution{V, F},
#                             z::ProbabilityDistribution{V, SampleList} = ProbabilityDistribution(V, SampleList)) where {V<:VariateType, F<:FactorNode}
#
#     thenode = currentGraph().nodes[x.params[:node_id]]
#
#     samples_ = [sample(prob, thenode.num_samples) for prob in thenode.q_memory]
#     samples = thenode.g.(samples_...)
#     weights = ones(thenode.num_samples)/thenode.num_samples
#     thenode.q_memory = thenode.q
#     return ProbabilityDistribution(Univariate, SampleList, s=samples, w=weights)
# end

@symmetrical function prod!(x::ProbabilityDistribution{Univariate, SetSampleList},
                            y::ProbabilityDistribution{Univariate, F},
                            z::ProbabilityDistribution{Univariate, SampleList} = ProbabilityDistribution(Univariate, SampleList)) where {F<:Gaussian}

    thenode = currentGraph().nodes[x.params[:node_id]]

    samples_ = [sample(prob, thenode.num_samples) for prob in thenode.q_memory]
    samples = thenode.g.(samples_...)
    weights = ones(thenode.num_samples)/thenode.num_samples
    thenode.q_memory = thenode.q
    return ProbabilityDistribution(Univariate, SampleList, s=samples, w=weights)
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, SetSampleList},
                            y::ProbabilityDistribution{Multivariate, F},
                            z::ProbabilityDistribution{Multivariate, SampleList} = ProbabilityDistribution(Multivariate, SampleList)) where {F<:Gaussian}

    thenode = currentGraph().nodes[x.params[:node_id]]

    samples_ = [sample(prob, thenode.num_samples) for prob in thenode.q_memory]
    samples = thenode.g.(samples_...)
    weights = ones(thenode.num_samples)/thenode.num_samples
    thenode.q_memory = thenode.q
    return ProbabilityDistribution(Multivariate, SampleList, s=samples, w=weights)
end

@symmetrical function prod!(x::ProbabilityDistribution{Multivariate, SetSampleList},
                            y::ProbabilityDistribution{MatrixVariate, Wishart},
                            z::ProbabilityDistribution{Multivariate, SampleList} = ProbabilityDistribution(Multivariate, SampleList))

    thenode = currentGraph().nodes[x.params[:node_id]]

    samples_ = [sample(prob, thenode.num_samples) for prob in thenode.q_memory]
    samples = thenode.g.(samples_...)
    weights = ones(thenode.num_samples)/thenode.num_samples
    thenode.q_memory = thenode.q
    return ProbabilityDistribution(Multivariate, SampleList, s=samples, w=weights)
end
