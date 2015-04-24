type QDistribution
    distribution::ProbabilityDistribution
    edges::Set{Edge} # Edges on which the distribution is defined
end

function calculateQDistribution!(node::Node, subgraph::Subgraph)
    # Calculate the approximate marginal for node from the perspective of subgraph,
    # and store the result in the scheme.q_distributions dictionary.

    q_distribution = current_algorithm.fields[:q_distributions][(node, subgraph)]
    if length(q_distribution.edges) == 1
        # Update for univariate q
        # When there is only one internal edge, the approximate marginal calculation reduces to the naive marginal update
        internal_edge = first(q_distribution.edges) # Extract element
        return calculateQDistribution!(q_distribution, internal_edge.tail.message.payload, internal_edge.head.message.payload)
    end

    # Update for multivariate q
    required_inputs = Array(Any, 0)
    for interface in node.interfaces # Iterate over all edges connected to node
        neighbouring_subgraph = current_algorithm.fields[:factorization].edge_to_subgraph[interface.edge]
        if neighbouring_subgraph == subgraph # edge is internal
            push!(required_inputs, interface.partner.message)
        else # edge is external
            haskey(current_algorithm.fields[:q_distributions], (node, neighbouring_subgraph)) || error("A required q-distribution for $(node.name) is not present. Please preset (vague) q-distributions.")
            push!(required_inputs, current_algorithm.fields[:q_distributions][(node, neighbouring_subgraph)])
        end
    end
    return calculateQDistribution!(q_distribution, node, required_inputs...)
end


############################
# Univariate approximate marginal calculations
############################

# GammaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    return gammaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# InverseGammaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# BetaDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::BetaDistribution, backward_dist::BetaDistribution)
    return betaMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# GaussianDistribution
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    return gaussianMarginalRule!(marg.distribution, forward_dist, backward_dist)
end

# Gaussian-students t combination
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    return gaussianStudentsMarginalRule!(marg.distribution, forward_dist, backward_dist)
end
calculateQDistribution!(marg::QDistribution, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(factor, graph, backward_dist, forward_dist)

# Gaussian-Delta combination
function calculateQDistribution!(marg::QDistribution, forward_dist::DeltaDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    return marg.distribution = convert(GaussianDistribution, forward_dist)
end
function calculateQDistribution!(marg::QDistribution, forward_dist::GaussianDistribution, backward_dist::DeltaDistribution)
    # Calculation for univariate approximate marginal
    return marg.distribution = convert(GaussianDistribution, backward_dist)
end


############################
# Joint approximate marginal calculations
############################

function calculateQDistribution!(q_dist::QDistribution,
                                 node::GaussianNode,
                                 gaus_msg::Message{GaussianDistribution},
                                 gam_msg::Message{GammaDistribution},
                                 y_dist::GaussianDistribution)
    # (Joint) marginal update function used for SVMP
    # Definitions available in derivations notebook

    marg = q_dist.distribution

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