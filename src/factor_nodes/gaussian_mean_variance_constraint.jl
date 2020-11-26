export GaussianMeanVarianceConstraint

"""
Description:

    Constraints the marginal of the connected variable to a Gaussian belief.
    
Interfaces:

    1. out

Construction:

    GaussianMeanVarianceConstraint(out; m=m, v=v, id=:my_node)
"""
mutable struct GaussianMeanVarianceConstraint <: SoftFactor # TODO: how to handle free energy evaluation for a point-mass constrained variable?
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    m::Union{Float64, AbstractVector}
    v::Union{Float64, AbstractMatrix}

    function GaussianMeanVarianceConstraint(out; m::Union{Float64, AbstractVector}, v::Union{Float64, AbstractMatrix},  id=ForneyLab.generateId(GaussianMeanVarianceConstraint))
        @ensureVariables(out)
        self = new(id, Vector{Interface}(undef, 1), Dict{Symbol,Interface}(), m, v)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end    

slug(::Type{GaussianMeanVarianceConstraint}) = "N_q"

# A breaker message is required if interface is partnered with a point-mass constraint
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::GaussianMeanVarianceConstraint) = true

function breakerParameters(interface::Interface, partner_interface::Interface, partner_node::GaussianMeanVarianceConstraint)
    return (Message{GaussianMeanVariance, variateType(partner_node.m)}, size(partner_node.m))
end

# Average energy functional
averageEnergy(::Type{GaussianMeanVarianceConstraint}, marg_out::ProbabilityDistribution) = differentialEntropy(marg_out)