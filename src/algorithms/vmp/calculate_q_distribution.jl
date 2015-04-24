function ensureQDistribution!(factor::Set{Edge}, scheme::DataAwareFactorGraph, assign_distribution::DataType)
    # Looks for a marginal in the q-distribution dictionary.
    # If no marginal is present, it sets and returns a vague distribution.
    # Otherwise, it returns the existing q distribution.
    try
        return qDistribution(factor)
    catch
        if assign_distribution <: ProbabilityDistribution
            return scheme.q_distributions[factor] = vague(assign_distribution) # Assign a q distribution
        else
            error("Cannot create a marginal of type $(assign_distribution) since a marginal should be <: ProbabilityDistribution")
        end
    end

end

# GammaDistribution
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GammaDistribution)
    return gammaMarginalRule!(marg, forward_dist, backward_dist)
end

# InverseGammaDistribution
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)
end

# BetaDistribution
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, BetaDistribution)
    return betaMarginalRule!(marg, forward_dist, backward_dist)
end

# GaussianDistribution
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GaussianDistribution)
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)
end

# Gaussian-students t combination
function calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    # Calculation for univariate approximate marginal
    marg = ensureQDistribution!(factor, scheme, GaussianDistribution)
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)
end
calculateQDistribution!(factor::Set{Edge}, scheme::InferenceScheme, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(factor, graph, backward_dist, forward_dist)

# Gaussian-Delta combination
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