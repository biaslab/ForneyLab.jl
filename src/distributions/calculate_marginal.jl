# Functions for calculating marginals on nodes and edges.

# Marginal calculations are in a separate file from the distribution definitions,
# because types that are needed for marginal calculations are defined after distribution definitions

# Marginal calculations on the edges are the same as the equality node rules,
# where the forward and backward messages are incoming and the marginal outcome is on the outgoing edge.

############################
# Edge marginal calculations
############################

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

# Gaussian-studens t combination
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    # Calculate the marginal from a forward/backward message pair.
    ensureMWParametrization!(forward_dist)
    (length(forward_dist.m) == 1 && length(forward_dist.W) == 1) || error("Equality node update for StudentsTDistribution and GaussianDistribution only supports univariate Gaussian distribution.")

    # Definitions available in derivations notebook
    l_a = backward_dist.W[1, 1] 
    mu_a = backward_dist.m[1]
    l_b = forward_dist.W[1, 1]
    mu_b = forward_dist.m[1]
    nu_term = 1 + (1/backward_dist.nu)

    return GaussianDistribution(m = (l_a*nu_term*mu_a + l_b*mu_b) / (l_a*nu_term + l_b), W = l_a*nu_term + l_b)
end
calculateMarginal(forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    marg = getOrCreateMarginal(edge, GaussianDistribution)
    ensureMWParametrization!(forward_dist)
    (length(forward_dist.m) == 1 && length(forward_dist.W) == 1) || error("Equality node update for StudentsTDistribution and GaussianDistribution only supports univariate Gaussian distribution.")

    # Definitions available in derivations notebook
    l_a = backward_dist.W[1, 1] 
    mu_a = backward_dist.m[1]
    l_b = forward_dist.W[1, 1]
    mu_b = forward_dist.m[1]
    nu_term = 1 + (1/backward_dist.nu)

    marg.m = [(l_a*nu_term*mu_a + l_b*mu_b) / (l_a*nu_term + l_b)]
    marg.W = reshape([l_a*nu_term + l_b], 1, 1)
    marg.xi = nothing
    marg.V = nothing
    return marg
end
calculateMarginal!(edge::Edge, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)

# Gaussian-Float64 combination
# Float64 can be seen as a Gaussian with zero variance (delta peak at the mean).
# Therefore a multiplication of a float with any Gaussian returns the float value.
calculateMarginal(forward_dist::Float64, ::GaussianDistribution) = deepcopy(forward_dist)
calculateMarginal(::GaussianDistribution, backward_dist::Float64) = deepcopy(backward_dist)
function calculateMarginal!(edge::Edge, forward_dist::Float64, ::GaussianDistribution)
    return edge.marginal = deepcopy(forward_dist)
end
function calculateMarginal!(edge::Edge, ::GaussianDistribution, backward_dist::Float64)
    return edge.marginal = deepcopy(backward_dist)
end


############################
# Node (joint) marginal calculations
############################

function calculateMarginal!(node::GaussianNode)
    # (Joint) marginal update function used for VMP
    # Definitions available in derivations notebook

    marg = getOrCreateMarginal(node, NormalGammaDistribution)

    mu_m = node.mean.partner.message.payload
    mu_gam = node.precision.partner.message.payload
    q_y = node.out.edge.marginal
    
    ensureMDefined!(mu_m)
    ensureMWParametrization!(q_y)

    (length(mu_m.m) == 1 && length(q_y.m) == 1)|| error("Update rule for NormalGammaDistribution marginal only supports univariate distributions.")

    marg.m = mu_m.m[1]
    marg.beta = 10000.0 # Big number
    marg.a = mu_gam.a + 0.5
    marg.b = (1.0/(2.0*q_y.W[1, 1])) + mu_gam.b

    return marg
end