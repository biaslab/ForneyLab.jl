function execute(subgraph::Subgraph, rf::RecognitionFactorization)
    # Execute internal schedule
    execute(subgraph.internal_schedule)

    # Update recognition distributions at external edges
    if ForneyLab.verbose
        println(" ")
        println("Marginals (node), result")
        println("--------------------------------------------")
    end
    for node in subgraph.nodes_connected_to_external_edges
        d = calculateRecognitionDistribution!(rf, node, subgraph)
        if ForneyLab.verbose
            println("$(node.id) : $(format(d))")
        end
    end
end

function execute(rf::RecognitionFactorization)
    for subgraph in rf.subgraphs
        execute(subgraph, rf)
    end
end

"""
Calculate the recognition distribution for the node-subgraph combination
"""
function calculateRecognitionDistribution!(rf::RecognitionFactorization, node::Node, subgraph::Subgraph)
    recognition_distribution = rf.node_subgraph_to_recognition_distribution[(node, subgraph)]
    internal_edges = rf.node_subgraph_to_internal_edges[(node, subgraph)]

    if length(internal_edges) == 1
        # Update for univariate q; when there is only one internal edge
        internal_edge = first(internal_edges)
        return prod!(internal_edge.tail.message.payload, internal_edge.head.message.payload, recognition_distribution)
    else
        # Update for multivariate q
        required_inputs = Array(Any, 0)
        for interface in node.interfaces # Iterate over all interfaces connected to node
            neighbouring_subgraph = rf.edge_to_subgraph[interface.edge]
            if neighbouring_subgraph == subgraph # edge is internal
                push!(required_inputs, interface.partner.message)
            else # edge is external
                push!(required_inputs, rf.node_subgraph_to_recognition_distribution[(node, neighbouring_subgraph)])
            end
        end
        return recognitionRule!(node, recognition_distribution, required_inputs...)
    end
end


############################
# Update rules for multivariate recognition distributions
############################

"""
NormalGamma update rule for variational message passing
"""
function recognitionRule!(  node::GaussianNode,
                            recognition_dist::NormalGamma,
                            msg_mean::Message{Gaussian},
                            msg_prec::Message{Gamma},
                            dist_y::Gaussian)

    dist_mean = msg_mean.payload
    dist_prec = msg_prec.payload
    ForneyLab.ensureParameters!(dist_mean, (:m,))
    ForneyLab.ensureParameters!(dist_y, (:W,))

    recognition_dist.m = dist_mean.m
    recognition_dist.beta = huge
    recognition_dist.a = dist_prec.a + 0.5
    recognition_dist.b = (1.0/(2.0*dist_y.W)) + dist_prec.b

    return recognition_dist
end