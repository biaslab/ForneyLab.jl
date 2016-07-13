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
        rdg = summaryDependencyGraph(sg, reverse_edges=true)

        # Find the interfaces whose messages are computed by the variational update rule.
        # All messages that are dependent on a variational message need to be part of the iterative schedule,
        # since the dependent messages might be altered during an iteration. Messages independent of a variational update
        # only need to be computer once, since they will never change during iterations.
        variational_interfaces = Interface[] # Interfaces for messages that need to be computed through the variational update rule
        partner_variational_interfaces = Interface[]
        for node in sg.nodes_connected_to_external_edges
            for interface in node.interfaces
                if interface.edge in sg.internal_edges
                    push!(variational_interfaces, interface) # The message on this interface must be computed by the variational update
                    push!(partner_variational_interfaces, interface.partner) # The message on this interface is not necessarily dependent on a variational update
                end
            end
        end
        influenced_by_variational_update = children(variational_interfaces, rdg)

        # TODO: this scheduler requires messages for the structured factorization that are not always required;
        # e.g. backward messages over the mean and precision edges present in a joint factorization when there
        # are no target dependencies.

        # TODO: when terminal nodes are also included in the nodes_connected_to_external_edges (which technically they don't belong to; except in case of a precision factor),
        # this might add messages dependent on prior nodes to the iterative schedule instead of the pre schedule (even though the prior node has a variational update,
        # the outbound message never changes; except in the case of a precision factor). 

        # Based on the variationally dependent messages, distinguish between the pre and iterative schedule
        pre_interface_list = Interface[]
        iterative_interface_list = Interface[]
        for interface in children(partner_variational_interfaces, dg)
            if interface in influenced_by_variational_update # Interface is influenced by a variational update, add to iterative schedule
                push!(iterative_interface_list, interface)
            else # Interface is not influenced by a variational update, only compute once in pre schedule
                push!(pre_interface_list, interface)
            end
        end
        iterative_interface_list = unique([iterative_interface_list; variational_interfaces]) # Add variational interfaces at the end because these may depend in internal messages in the case of a structured factorization.

        # Propagate messages towards targets (such as wraps and write buffers)
        post_interface_list = Interface[]
        for interface in targets
            if interface.edge in sg.internal_edges # target is the responsibility of the subgraph sg
                post_interface_list = [post_interface_list; children(interface, dg, breakers=union(Set(pre_interface_list), Set(iterative_interface_list)))]
            end
        end
        
        # Covert interface lists to schedule; assign variational update rules to variational interfaces
        sg.internal_pre_schedule = convertToVariationalSchedule(pre_interface_list, variational_interfaces)
        sg.internal_iterative_schedule = convertToVariationalSchedule(iterative_interface_list, variational_interfaces)
        sg.internal_post_schedule = convertToVariationalSchedule(post_interface_list, variational_interfaces)

        # Infer types on pre schedule
        skip_update_interfaces = Interface[] # Interfaces for which the variational update need not be computed
        for entry in sg.internal_pre_schedule
            inferTypes!(entry, message_types, recognition_factorization)

            # Find interfaces of superfluous variational updates
            outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
            if (outbound_interface in partner_variational_interfaces) && (entry.outbound_type <: AbstractDelta)
                # The variational message will be combined with a delta message to yield the recognition distribution.
                # In this case the recognition distribution will become equal to the delta message,
                # so it is superfluous to compute the variational message.
                push!(skip_update_interfaces, outbound_interface.partner)
            end
        end

        # Infer types on iterative schedule
        for entry in sg.internal_iterative_schedule
            # Find schedule entries of superfluous variational updates
            outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
            # TODO: optionally turn skipping off
            if outbound_interface in skip_update_interfaces
                deleteat!(sg.internal_iterative_schedule, findfirst(sg.internal_iterative_schedule, entry)) # Remove entry from schedule
            else
                inferTypes!(entry, message_types, recognition_factorization)
            end
        end

        # Infer types on post schedule
        for entry in sg.internal_post_schedule
            inferTypes!(entry, message_types, recognition_factorization)
        end
    end

    algo = VariationalBayes(graph, n_iterations, recognition_factorization)

    return algo
end

function convertToVariationalSchedule(interface_list::Vector{Interface}, variational_interfaces::Vector{Interface})
    schedule = ScheduleEntry[]
    for interface in interface_list
        if interface in variational_interfaces
            push!(schedule, convert(ScheduleEntry{VariationalRule}, interface))
        else
            push!(schedule, convert(ScheduleEntry{SumProductRule}, interface))
        end
    end

    return schedule
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
        schedules = (factor.internal_pre_schedule, factor.internal_iterative_schedule, factor.internal_post_schedule)

        # Populate the subgraph with vague messages of the correct types
        for schedule in schedules
            for entry in schedule
                ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
            end
        end

        # Compile the schedules (define schedule_entry.execute)
        for schedule in schedules
            compile!(schedule, algo)
        end
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
    rf = algo.recognition_factorization
    resetRecognitionDistributions!(rf)

    # Execute internal pre schedules
    for subgraph in rf.subgraphs
        isempty(subgraph.internal_pre_schedule) || execute(subgraph.internal_pre_schedule)
    end

    # Execute iterative schedules
    for iteration = 1:algo.n_iterations
        for subgraph in rf.subgraphs
            # Execute internal schedule
            isempty(subgraph.internal_iterative_schedule) || execute(subgraph.internal_iterative_schedule)

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
    end

    # Execute internal post schedules
    for subgraph in rf.subgraphs
        isempty(subgraph.internal_post_schedule) || execute(subgraph.internal_post_schedule)
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
        (internal_edge.tail.message == nothing) ? tail_dist = nothing : tail_dist = internal_edge.tail.message.payload
        (internal_edge.head.message == nothing) ? head_dist = nothing : head_dist = internal_edge.head.message.payload
        return prod!(tail_dist, head_dist, recognition_distribution)
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