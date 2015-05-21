export calculateMarginal, calculateMarginal!, getMarginalType

# Functions for calculating marginals on nodes and edges.

# Marginal calculations are in a separate file from the distribution definitions,
# because types that are needed for marginal calculations are defined after distribution definitions

# Marginal calculations on the edges are the same as the equality node rules,
# where the forward and backward messages are incoming and the marginal outcome is on the outgoing edge.

############################
# Edge marginal calculations
############################

function calculateMarginal(edge::Edge)
    # Calculates the marginal without writing back to the edge
    @assert(edge.tail.message != nothing, "Edge ($(edge.tail.node.name) --> $(edge.head.node.name)) should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge ($(edge.tail.node.name) --> $(edge.head.node.name)) should hold a backward message.")
    return calculateMarginal(edge.tail.message.payload, edge.head.message.payload)
end
function calculateMarginal!(edge::Edge)
    # Calculates and writes the marginal on edge
    @assert(edge.tail.message != nothing, "Edge ($(edge.tail.node.name) --> $(edge.head.node.name)) should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge ($(edge.tail.node.name) --> $(edge.head.node.name)) should hold a backward message.")
    calculateMarginal!(edge, edge.tail.message.payload, edge.head.message.payload)
    return edge.marginal
end

# GammaDistribution
function calculateMarginal(forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    marg = GammaDistribution() # Do not overwrite an existing distribution
    return equalityGammaRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    marg = ensureMarginal!(edge, GammaDistribution)
    return equalityGammaRule!(marg, forward_dist, backward_dist)
end

# InverseGammaDistribution   
function calculateMarginal(forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    marg = InverseGammaDistribution() # Do not overwrite an existing distribution
    return equalityInverseGammaRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    marg = ensureMarginal!(edge, InverseGammaDistribution)
    return equalityInverseGammaRule!(marg, forward_dist, backward_dist)
end

# BetaDistribution
function calculateMarginal(forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    marg = BetaDistribution() # Do not overwrite an existing distribution
    return equalityBetaRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    marg = ensureMarginal!(edge, BetaDistribution)
    return equalityBetaRule!(marg, forward_dist, backward_dist)
end

# GaussianDistribution
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    marg = GaussianDistribution() # Do not overwrite an existing distribution
    return equalityGaussianRule!(marg, forward_dist, backward_dist)  
end
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    marg = ensureMarginal!(edge, GaussianDistribution)
    return equalityGaussianRule!(marg, forward_dist, backward_dist)
end

# DeltaDistribution
function calculateMarginal(forward_dist::DeltaDistribution, backward_dist::DeltaDistribution)
    marg = DeltaDistribution() # Do not overwrite an existing distribution
    return equalityDeltaRule!(marg, forward_dist, backward_dist)  
end
function calculateMarginal!(edge::Edge, forward_dist::DeltaDistribution, backward_dist::DeltaDistribution)
    marg = ensureMarginal!(edge, DeltaDistribution)
    return equalityDeltaRule!(marg, forward_dist, backward_dist)
end

# Gaussian-students t combination
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    marg = GaussianDistribution() # Do not overwrite an existing distribution
    return equalityGaussianStudentsTRule!(marg, forward_dist, backward_dist)    
end
calculateMarginal(forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    marg = ensureMarginal!(edge, GaussianDistribution)
    return equalityGaussianStudentsTRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)

# Gaussian-Delta combination
# A multiplication of a delta distribution with any Gaussian returns the delta.
function calculateMarginal(forward_dist::DeltaDistribution, backward_dist::GaussianDistribution)
    marg = DeltaDistribution()
    return equalityGaussianDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal(forward_dist::GaussianDistribution, backward_dist::DeltaDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::DeltaDistribution, backward_dist::GaussianDistribution)
    marg = ensureMarginal!(edge, DeltaDistribution)
    return equalityGaussianDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)

# Gamma-Delta combination
# A multiplication of a delta distribution (with a peak >= 0) with any gamma returns the delta.
function calculateMarginal(forward_dist::DeltaDistribution, backward_dist::GammaDistribution)
    marg = DeltaDistribution()
    return equalityGammaDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal(forward_dist::GammaDistribution, backward_dist::DeltaDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::DeltaDistribution, backward_dist::GammaDistribution)
    marg = ensureMarginal!(edge, DeltaDistribution)
    return equalityGammaDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::GammaDistribution, backward_dist::DeltaDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)

# InverseGamma-Delta combination
# A multiplication of a delta distribution (with a peak >= 0) with any inverse gamma returns the delta.
function calculateMarginal(forward_dist::DeltaDistribution, backward_dist::InverseGammaDistribution)
    marg = DeltaDistribution()
    return equalityInverseGammaDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal(forward_dist::InverseGammaDistribution, backward_dist::DeltaDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::DeltaDistribution, backward_dist::InverseGammaDistribution)
    marg = ensureMarginal!(edge, DeltaDistribution)
    return equalityInverseGammaDeltaRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::InverseGammaDistribution, backward_dist::DeltaDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)

########################################
# Lookup table for joint marginals
########################################
getMarginalType(::Type{GaussianDistribution},   ::Type{GammaDistribution})      = NormalGammaDistribution
getMarginalType(::Type{GammaDistribution},      ::Type{GaussianDistribution})   = NormalGammaDistribution
