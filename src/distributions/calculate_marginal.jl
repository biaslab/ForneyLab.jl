export calculateMarginal, calculateMarginal!, calculateQDistribution!

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
function gammaMarginalRule!(marg::GammaDistribution, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg.a = forward_dist.a+backward_dist.a-1.0
    marg.b = forward_dist.b+backward_dist.b
    return marg
end    
function calculateMarginal(forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    marg = GammaDistribution() # Do not overwrite an existing distribution
    return gammaMarginalRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    marg = ensureMarginal!(edge, GammaDistribution)
    return gammaMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GammaDistribution)
    return gammaMarginalRule!(marg, forward_dist, backward_dist)
end

# InverseGammaDistribution
function inverseGammaMarginalRule!(marg::InverseGammaDistribution, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg.a = forward_dist.a+backward_dist.a+1.0
    marg.b = forward_dist.b+backward_dist.b
    return marg
end    
function calculateMarginal(forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    marg = InverseGammaDistribution() # Do not overwrite an existing distribution
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    marg = ensureMarginal!(edge, InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)
end

# BetaDistribution
function betaMarginalRule!(marg::BetaDistribution, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    marg.a = forward_dist.a+backward_dist.a-1.0
    marg.b = forward_dist.b+backward_dist.b-1.0
    return marg
end    
function calculateMarginal(forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    marg = BetaDistribution() # Do not overwrite an existing distribution
    return betaMarginalRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    marg = ensureMarginal!(edge, BetaDistribution)
    return betaMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, BetaDistribution)
    return betaMarginalRule!(marg, forward_dist, backward_dist)
end

# GaussianDistribution
function gaussianMarginalRule!(marg::GaussianDistribution, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    ensureXiWParametrization!(forward_dist)
    ensureXiWParametrization!(backward_dist)
    marg.xi = forward_dist.xi+backward_dist.xi
    marg.W = forward_dist.W+backward_dist.W
    marg.V = nothing
    marg.m = nothing
    return marg
end    
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    marg = GaussianDistribution() # Do not overwrite an existing distribution
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)    
end
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    marg = ensureMarginal!(edge, GaussianDistribution)
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GaussianDistribution)
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)
end

# Gaussian-students t combination
function gaussianStudentsMarginalRule!(marg::GaussianDistribution, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using the EqualityNode update rules; same for the functions below
    ensureMWParametrization!(forward_dist)
    (length(forward_dist.m) == 1 && length(forward_dist.W) == 1) || error("Marginal update for StudentsTDistribution and GaussianDistribution only supports univariate Gaussian distribution.")

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
function calculateMarginal(forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    marg = GaussianDistribution() # Do not overwrite an existing distribution
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)    
end
calculateMarginal(forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal(backward_dist, forward_dist)
function calculateMarginal!(edge::Edge, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    marg = ensureMarginal!(edge, GaussianDistribution)
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GaussianDistribution)
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)
end
calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(factor, graph, backward_dist, forward_dist)

# Gaussian-Delta combination
# A multiplication of a delta distribution with any Gaussian returns the delta.
calculateMarginal(forward_dist::DeltaDistribution, ::GaussianDistribution) = deepcopy(forward_dist)
calculateMarginal(::GaussianDistribution, backward_dist::DeltaDistribution) = deepcopy(backward_dist)
function calculateMarginal!(edge::Edge, forward_dist::DeltaDistribution, ::GaussianDistribution)
    return edge.marginal = convert(GaussianDistribution, forward_dist)
end
function calculateMarginal!(edge::Edge, ::GaussianDistribution, backward_dist::DeltaDistribution)
    return edge.marginal = convert(GaussianDistribution, backward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::DeltaDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    return scheme.q_distributions[factor] = convert(GaussianDistribution, forward_dist)
end
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution)
    # Calculation for univariate approximate marginal
    return scheme.q_distributions[factor] = convert(GaussianDistribution, backward_dist)
end


############################
# Joint approximate marginal calculations
############################

function calculateQDistribution!(node::Node, factor::Set{Edge}, scheme::InferenceScheme)
    # Calculate the approximate marginal for node from the perspective of subgraph,
    # and store the result in the scheme.q_distributions dictionary.

    if length(factor) == 1
        # Update for univariate q
        # When there is only one internal edge, the approximate marginal calculation reduces to the naive marginal update
        internal_edge = first(factor) # Extract element
        calculateQDistribution!(factor, scheme, internal_edge.tail.message.payload, internal_edge.head.message.payload)
        return qDistribution(scheme, factor)
    end

    # Update for multivariate q
    required_inputs = Array(Any, 0)
    local_factors = qFactors(scheme, node)
    for j = 1:length(node.interfaces) # Iterate over all edges connected to node
        if factor == local_factors[j] # edge is internal
            push!(required_inputs, node.interfaces[j].partner.message)
        else # edge is external
            haskey(scheme.q_distributions, local_factors[j]) || error("A required approximate marginal for $(node.name) is not preset. Please preset an (vague) marginal.")
            push!(required_inputs, qDistribution(scheme, local_factors[j]))
        end
    end
    calculateQDistribution!(factor, scheme, required_inputs...)
    return qDistribution(scheme, factor)
end

function calculateQDistribution!(factor::Set{Edge},
                            scheme::InferenceScheme,
                            gaus_msg::Message{GaussianDistribution},
                            gam_msg::Message{GammaDistribution},
                            y_dist::GaussianDistribution)
    # (Joint) marginal update function used for SVMP
    # Definitions available in derivations notebook

    marg = ensureQDistribution!(factor, scheme, NormalGammaDistribution)

    mu_m = gaus_msg.payload
    mu_gam = gam_msg.payload
    
    ensureMDefined!(mu_m)
    ensureMWParametrization!(y_dist)

    (length(mu_m.m) == 1 && length(y_dist.m) == 1)|| error("Update rule for NormalGammaDistribution marginal only supports univariate distributions.")

    marg.m = mu_m.m[1]
    marg.beta = huge()
    marg.a = mu_gam.a + 0.5
    marg.b = (1.0/(2.0*y_dist.W[1, 1])) + mu_gam.b

    return marg
end


########################################
# Lookup table for joint marginals
########################################
getMarginalType(::Type{GaussianDistribution},   ::Type{GammaDistribution})      = NormalGammaDistribution
getMarginalType(::Type{GammaDistribution},      ::Type{GaussianDistribution})   = NormalGammaDistribution
