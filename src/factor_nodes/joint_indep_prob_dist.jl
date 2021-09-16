export JointIndependentProbDist

"""
Description:

    JointIndependentProbDist is joint probability dist. of independent random variables.
    Useful to calculate differential entropy around CVI node with more than one inputs.
    Around CVI node, we assume a fully factorized recognition distribution.

"""
mutable struct JointIndependentProbDist <: SoftFactor
    params::Dict
end

slug(::Type{JointIndependentProbDist}) = "jipd"

format(dist::ProbabilityDistribution{Multivariate, JointIndependentProbDist}) = "$(dist.params)"

# Distribution constructors
ProbabilityDistribution(::Type{Multivariate}, ::Type{JointIndependentProbDist}; kwargs...) = ProbabilityDistribution{Multivariate, SetProbDist}(Dict{Symbol,Any}(kwargs))
ProbabilityDistribution(::Type{JointIndependentProbDist}; kwargs...) = ProbabilityDistribution{Multivariate, JointIndependentProbDist}(Dict{Symbol,Any}(kwargs))

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Multivariate, JointIndependentProbDist})
    res = 0
    for marg in dist.params[:marginals]
        res += differentialEntropy(marg)
    end
    return res
end