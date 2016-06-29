import Base.show
export VariationalBayes

include("subgraph.jl")
include("recognition_factorization.jl")

"""
Variational message passing algorithm.

Usage:

    VariationalBayes(recognition_distribution_types::Dict, graph::Graph; n_iterations, message_types)
"""
type VariationalBayes <: InferenceAlgorithm
    graph::FactorGraph
    n_iterations::Int64
    recognition_factorization::RecognitionFactorization
end

function show(io::IO, algo::VariationalBayes)
    println("VariationalBayes inference algorithm")
    println("    number of subgraphs: $(length(algo.recognition_factorization.subgraphs))")
    println("    number of iterations: $(algo.n_iterations)")
end

############################################
# VariationalBayes algorithm constructors
############################################

"""
Generate a VariationalBayes algorithm
"""
function VariationalBayes(  graph::FactorGraph=currentGraph(),
                            recognition_factorization::RecognitionFactorization=currentRecognitionFactorization();
                            kwargs...)

    VariationalBayes(interfacesFacingWrapsOrBuffers(graph), graph, recognition_factorization; kwargs...)
end

function VariationalBayes(  targets::Vector{Interface},
                            graph::FactorGraph=currentGraph(),
                            recognition_factorization::RecognitionFactorization=currentRecognitionFactorization();
                            n_iterations::Int64=50,
                            message_types::Dict{Interface,DataType}=Dict{Interface,DataType}())
    
    verifyProper(recognition_factorization)

    # Generate the internal schedule for each subgraph sg
    for sg in recognition_factorization.subgraphs
        dg = summaryDependencyGraph(sg) # Dependency graph of subgraph internal edges
        
        # TODO: this scheduler requires messages for the structured factorization that are not always required;
        # e.g. backward messages over the mean and precision edges present in a joint factorization when there
        # are no target dependencies.

        # Update the messages along the internal edges of sg
        variational_interfaces = Interface[]
        internal_interface_list = Interface[]
        for g_node in sg.nodes_connected_to_external_edges # For each node g connected to external edges
            for outbound_interface in g_node.interfaces # Iterate over the outbound interfaces of the external node (messages on these interfaces will need to be computer through the variational rule)
                if outbound_interface.edge in sg.internal_edges # If the edge is internal, compute the forward and backward message
                    push!(variational_interfaces, outbound_interface) # Store sites for variational updates
                    internal_interface_list = [internal_interface_list; children([outbound_interface, outbound_interface.partner], dg)] # Build interface array for schedule
                end
            end
        end
        
        # Propagate messages towards targets (such as wraps and write buffers)
        target_interface_list = Interface[]
        for interface in targets
            if interface.edge in sg.internal_edges # target is the responsibility of the subgraph sg
                target_interface_list = [target_interface_list; children(interface, dg)]
            end
        end
        
        # Covert interface lists to schedule
        internal_schedule = ScheduleEntry[]
        for interface in unique([internal_interface_list; target_interface_list])
            if interface in variational_interfaces
                push!(internal_schedule, convert(ScheduleEntry{VariationalRule}, interface))
            else
                push!(internal_schedule, convert(ScheduleEntry{SumProductRule}, interface))
            end
        end

        # Assign internal schedule
        sg.internal_schedule = internal_schedule

        # Execute type inference
        for entry in sg.internal_schedule
            inferTypes!(entry, message_types, recognition_factorization)
        end
    end

    algo = VariationalBayes(graph, n_iterations, recognition_factorization)

    return algo    
end


############################################
# Type inference and preparation
############################################

inferTypes!(entry::ScheduleEntry{SumProductRule}, inferred_outbound_types::Dict{Interface, DataType}, ::RecognitionFactorization) = inferTypes!(entry, inferred_outbound_types) # Default to sum-product rule

function inferTypes!(entry::ScheduleEntry{VariationalRule}, inferred_outbound_types::Dict{Interface, DataType}, rf::RecognitionFactorization)
    outbound_interface = entry.node.interfaces[entry.outbound_interface_id]

    # Collect inbound types
    entry.inbound_types = DataType[]
    for (id, interface) in enumerate(entry.node.interfaces)
        sg = rf.edge_to_subgraph[interface.edge]
        if id == entry.outbound_interface_id
            push!(entry.inbound_types, Void)
        elseif is(sg, rf.edge_to_subgraph[outbound_interface.edge])
            # Both edges are in the same subgraph (structured factorization), require message
            push!(entry.inbound_types, Message{inferred_outbound_types[interface.partner]})
        else
            # A subgraph border is crossed, require recognition distribution
            push!(entry.inbound_types, typeof(rf.node_subgraph_to_recognition_distribution[(entry.node, sg)]))
        end
    end

    # Infer outbound type
    if outbound_interface in keys(inferred_outbound_types)
        setOutboundType!(entry, inferred_outbound_types[outbound_interface])
    end
    inferOutboundType!(entry) # If entry.outbound_type is already set, this will validate that there is a suitable rule available
    inferred_outbound_types[outbound_interface] = entry.outbound_type

    return entry
end

function prepare!(algo::VariationalBayes)
    for factor in algo.recognition_factorization.subgraphs
        schedule = factor.internal_schedule

        # Populate the subgraph with vague messages of the correct types
        for entry in schedule
            ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
        end

        # Compile the schedule (define schedule_entry.execute)
        compile!(schedule, algo)
    end

    return algo.graph.prepared_algorithm = algo
end

function compile!(entry::ScheduleEntry{VariationalRule}, algo::VariationalBayes)
    rf = algo.recognition_factorization
    outbound_interface = entry.node.interfaces[entry.outbound_interface_id]

    inbounds = []
    for (id, interface) in enumerate(entry.node.interfaces)
        sg = rf.edge_to_subgraph[interface.edge]
        if id == entry.outbound_interface_id
            push!(inbounds, rf.node_subgraph_to_recognition_distribution[(entry.node, sg)])
        elseif is(sg, rf.edge_to_subgraph[outbound_interface.edge])
            # Both edges are in the same subgraph (structured factorization), require message
            push!(inbounds, interface.partner.message)
        else
            # A subgraph border is crossed, require recognition distribution
            push!(inbounds, rf.node_subgraph_to_recognition_distribution[(entry.node, sg)])
        end
    end

    return buildExecute!(entry, inbounds)
end

function execute(algo::VariationalBayes)
    # Make sure algo is prepared
    (algo.graph.prepared_algorithm == algo) || prepare!(algo)

    # Reset recognition distributions to vague before next step
    resetRecognitionDistributions!(algo.recognition_factorization)

    for iteration = 1:algo.n_iterations
        execute(algo.recognition_factorization)
    end
end

function execute(rf::RecognitionFactorization)
    for subgraph in rf.subgraphs
        execute(subgraph, rf)
    end
end

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

"""
Calculate the recognition distribution for the node-subgraph combination
"""
function calculateRecognitionDistribution!(rf::RecognitionFactorization, node::Node, subgraph::Subgraph)
    recognition_distribution = rf.node_subgraph_to_recognition_distribution[(node, subgraph)]
    internal_edges = rf.node_subgraph_to_internal_edges[(node, subgraph)]

    if length(internal_edges) == 1
        # Update for univariate recognition distribution
        internal_edge = first(internal_edges)
        return prod!(internal_edge.tail.message.payload, internal_edge.head.message.payload, recognition_distribution)
    else
        # Update for multivariate recognition distribution
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