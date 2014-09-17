export calculateMarginal, calculateMarginal!

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
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
    return calculateMarginal(edge.tail.message.payload, edge.head.message.payload)
end
function calculateMarginal!(edge::Edge)
    # Calculates and writes the marginal on edge
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
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
    marg = getOrCreateMarginal(edge, GammaDistribution)
    return gammaMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::GammaDistribution, backward_dist::GammaDistribution)
    # Calculation for univariate approximate marginal
    marg = getOrCreateMarginal(node, subgraph, graph, GammaDistribution)
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
    marg = getOrCreateMarginal(edge, InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::InverseGammaDistribution, backward_dist::InverseGammaDistribution)
    # Calculation for univariate approximate marginal
    marg = getOrCreateMarginal(node, subgraph, graph, InverseGammaDistribution)
    return inverseGammaMarginalRule!(marg, forward_dist, backward_dist)
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
    marg = getOrCreateMarginal(edge, GaussianDistribution)
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)
end
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::GaussianDistribution, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    marg = getOrCreateMarginal(node, subgraph, graph, GaussianDistribution)
    return gaussianMarginalRule!(marg, forward_dist, backward_dist)
end

# Gaussian-studens t combination
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
    marg = getOrCreateMarginal(edge, GaussianDistribution)
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(edge::Edge, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(edge, backward_dist, forward_dist)
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::GaussianDistribution, backward_dist::StudentsTDistribution)
    # Calculation for univariate approximate marginal
    marg = getOrCreateMarginal(node, subgraph, graph, GaussianDistribution)
    return gaussianStudentsMarginalRule!(marg, forward_dist, backward_dist)
end
calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::StudentsTDistribution, backward_dist::GaussianDistribution) = calculateMarginal!(node, subgraph, graph, backward_dist, forward_dist)

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
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::Float64, backward_dist::GaussianDistribution)
    # Calculation for univariate approximate marginal
    return graph.approximate_marginals[(node, subgraph)] = deepcopy(forward_dist)
end
function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, forward_dist::GaussianDistribution, backward_dist::Float64)
    # Calculation for univariate approximate marginal
    return graph.approximate_marginals[(node, subgraph)] = deepcopy(backward_dist)
end


############################
# Joint approximate marginal calculations
############################

function calculateMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph=getCurrentGraph())
    # Calculate the approximate marginal for node from the perspective of subgraph,
    # and store the result in the graph.approximate_marginals dictionary.

    # Gather internal and external messages/qs
    required_inputs = Array(Union(Message, MessagePayload), 0)
    internal_edge_list = Array(Edge, 0)
    for interface in node.interfaces # In the order of the node's interfaces
        neighbouring_subgraph = graph.edge_to_subgraph[interface.edge]
        if neighbouring_subgraph == subgraph # Edge is internal
            push!(required_inputs, interface.partner.message)
            push!(internal_edge_list, interface.edge)
        else # Edge is external
            haskey(graph.approximate_marginals, (node, neighbouring_subgraph)) || error("A required approximate marginal for $(node.name) is not preset. Please preset an (uninformative) marginal.")
            push!(required_inputs, graph.approximate_marginals[(node, neighbouring_subgraph)])
        end
    end

    if length(internal_edge_list) == 1
        # Update for univariate q
        # When there is only one internal edge, the approximate marginal calculation reduces to the naive marginal update
        # The internal_edge_list contains the edge that requires the update
        internal_edge = internal_edge_list[1]
        calculateMarginal!(node, subgraph, graph, internal_edge.tail.message.payload, internal_edge.head.message.payload) 
    else
        # Update for multivariate q
        # This update involves the node function 
        calculateMarginal!(node, subgraph, graph, required_inputs...)
    end

    return graph.approximate_marginals[(node, subgraph)]
end

function calculateMarginal!(node::GaussianNode,
                            subgraph::Subgraph,
                            graph::FactorGraph,
                            gaus_msg::Message{GaussianDistribution},
                            gam_msg::Message{GammaDistribution},
                            y_dist::GaussianDistribution)
    # (Joint) marginal update function used for SVMP
    # Definitions available in derivations notebook

    marg = getOrCreateMarginal(node, subgraph, graph, NormalGammaDistribution)

    mu_m = gaus_msg.payload
    mu_gam = gam_msg.payload
    
    ensureMDefined!(mu_m)
    ensureMWParametrization!(y_dist)

    (length(mu_m.m) == 1 && length(y_dist.m) == 1)|| error("Update rule for NormalGammaDistribution marginal only supports univariate distributions.")

    marg.m = mu_m.m[1]
    marg.beta = 10000.0 # Big number
    marg.a = mu_gam.a + 0.5
    marg.b = (1.0/(2.0*y_dist.W[1, 1])) + mu_gam.b

    return marg
end


############################
# Lookup for marginal type combinations
############################

# Univariate
getMarginalType(::Type{GaussianDistribution}, ::Type{GaussianDistribution}) = GaussianDistribution
getMarginalType(::Type{GammaDistribution}, ::Type{GammaDistribution}) = GammaDistribution
getMarginalType(::Type{InverseGammaDistribution}, ::Type{InverseGammaDistribution}) = InverseGammaDistribution
getMarginalType(::Type{GaussianDistribution}, ::Type{StudentsTDistribution}) = GaussianDistribution
getMarginalType(::Type{StudentsTDistribution}, ::Type{GaussianDistribution}) = GaussianDistribution
getMarginalType(::Type{GaussianDistribution}, ::Type{Float64}) = Float64
getMarginalType(::Type{Float64}, ::Type{GaussianDistribution}) = Float64

# Multivariate
getMarginalType(::Type{GaussianNode}, ::Type{GaussianDistribution}, ::Type{GammaDistribution}) = NormalGammaDistribution
getMarginalType(::Type{GaussianNode}, ::Type{GammaDistribution}, ::Type{GaussianDistribution}) = NormalGammaDistribution
