# Helper functions for distributions

#######################
# Marginal calculations
#######################

# GammaDistribution
function calculateMarginal(forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    return GammaDistribution(a = forward_dist.a+backward_dist.a-1.0, b = forward_dist.b+backward_dist.b)    
end
function calculateMarginal!(edge::Edge, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg = getOrCreateMarginal(edge, GammaDistribution)
    marg.a = forward_dist.a+backward_dist.a-1.0
    marg.b = forward_dist.b+backward_dist.b
    return marg
end

# InverseGammaDistribution
function calculateMarginal(forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    return InverseGammaDistribution(a = forward_dist.a+backward_dist.a+1.0, b = forward_dist.b+backward_dist.b)    
end
function calculateMarginal!(edge::Edge, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg = getOrCreateMarginal(edge, InverseGammaDistribution)
    marg.a = forward_dist.a+backward_dist.a+1.0
    marg.b = forward_dist.b+backward_dist.b
    return marg
end

# GaussianDistribution
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    ensureXiWParametrization!(forward_dist)
    ensureXiWParametrization!(backward_dist)
    return GaussianDistribution(xi = forward_dist.xi+backward_dist.xi, W = forward_dist.W+backward_dist.W)
end
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg = getOrCreateMarginal(edge, GaussianDistribution)
    ensureXiWParametrization!(forward_dist)
    ensureXiWParametrization!(backward_dist)
    marg.xi = forward_dist.xi+backward_dist.xi
    marg.W = forward_dist.W+backward_dist.W
    marg.V = nothing
    marg.m = nothing
    return marg
end