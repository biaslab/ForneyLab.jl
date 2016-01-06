import Base.show
export VariationalBayes

include("subgraph.jl")
include("q_distribution.jl")
include("q_factorization.jl")
include("scheduler.jl")
include("message_passing.jl")

type VariationalBayes <: InferenceAlgorithm
    execute::Function
    factorization::QFactorization
    q_distributions::Dict{Tuple{Node,Subgraph},QDistribution}
    n_iterations::Int64
end


############################################
# VariationalBayes algorithm constructors
############################################

function VariationalBayes(graph::FactorGraph=currentGraph(); n_iterations::Int64=50)
    # Generates a VariationalBayes algorithm that propagates messages to all write buffers and wraps.
    # Uses a mean field factorization and autoscheduler
    factorization = factorize(graph) # Mean field factorization
    generateVariationalBayesSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = initializeVagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        resetQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.factorization, algorithm.q_distributions) # For all subgraphs, execute internal and external schedules
        end
    end

    algo = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    inferDistributionTypes!(algo)

    return algo
end

function VariationalBayes(graph::FactorGraph, cluster_edges...; n_iterations::Int64=50)
    # Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on graph
    # Uses a structured factorization with elements in cluster_edges defining separate subgraphs
    factorization = QFactorization(graph)
    for cluster in cluster_edges
        factorization = factorize!(cluster, factorization) # Structured factorization with elements in cluster_edges defining separate subgraphs
    end
    generateVariationalBayesSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = initializeVagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        resetQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.factorization, algorithm.q_distributions) # For all subgraphs, execute internal and external schedules
        end
    end

    algo = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    inferDistributionTypes!(algo)

    return algo
end
VariationalBayes(cluster_edges...; n_iterations::Int64=50) = VariationalBayes(currentGraph(), cluster_edges...; n_iterations=n_iterations)


############################################
# Type inference and preparation 
############################################

function inferDistributionTypes!(algo::VariationalBayes)
    # Infer the payload types for all messages in the internal schedules

    for factor in algo.factorization.factors
        # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
        schedule = factor.internal_schedule
        schedule_entries = Dict{Interface, ScheduleEntry}()

        for entry in schedule
            # Generate array of inbound types
            node = entry.node
            outbound_interface_id = entry.outbound_interface_id
            outbound_interface = node.interfaces[outbound_interface_id]

            collectInboundTypes!(entry, node, outbound_interface_id, schedule_entries, algo) # VariationalBayes algorithm specific collection of inbound types
            inferOutboundType!(entry, node, outbound_interface_id, entry.inbound_types, [sumProduct!, vmp!]) # The VariationalBayes algorithm allows access to sumProduct! and vmp! update rules

            schedule_entries[outbound_interface] = entry # Assign schedule entry to lookup dictionary
        end
    end

    return algo
end

function collectInboundTypes!(entry::ScheduleEntry, node::Node, outbound_interface_id::Int64, schedule_entries::Dict{Interface, ScheduleEntry}, algo::VariationalBayes)
    entry.inbound_types = []
    outbound_interface = node.interfaces[outbound_interface_id]

    # Collect references to all required inbound messages for executing message computation rule
    for (id, interface) in enumerate(node.interfaces)
        if id == outbound_interface_id
            push!(entry.inbound_types, Void) # This interface is outbound, push Void
        else
            # Should we require the inbound message or marginal?
            if is(algo.factorization.edge_to_subgraph[interface.edge], algo.factorization.edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = algo.factorization.edge_to_subgraph[interface.edge]
                push!(entry.inbound_types, typeof(algo.q_distributions[(node, sg)].distribution))
            end
        end
    end

    return entry
end

function prepare!(algo::VariationalBayes)
    for factor in algo.factorization.factors
        schedule = factor.internal_schedule

        # Populate the subgraph with vague messages of the correct types
        for entry in schedule
            ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
        end

        # Compile the schedule (define schedule_entry.execute)
        compile!(schedule, algo)
    end

    return algo
end

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol("ForneyLab.vmp!")}}, algo::VariationalBayes)
    # Generate schedule_entry.execute for schedule entry with vmp update rule

    # Collect references to all required inbound messages for executing message computation rule
    node = schedule_entry.node
    outbound_interface_id = schedule_entry.outbound_interface_id
    outbound_interface = node.interfaces[outbound_interface_id]

    rule_arguments = []
    # Add inbound messages to rule_arguments
    for (id, interface) in enumerate(node.interfaces)
        if id == outbound_interface_id
            push!(rule_arguments, nothing) # This interface is outbound, push "nothing"
        else
            # Should we require the inbound message or marginal?
            if is(algo.factorization.edge_to_subgraph[interface.edge], algo.factorization.edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                push!(rule_arguments, interface.partner.message)
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = algo.factorization.edge_to_subgraph[interface.edge]
                push!(rule_arguments, algo.q_distributions[(node, sg)].distribution)
            end
        end
    end
    # Add outbound distribution to rule_arguments
    push!(rule_arguments, node.interfaces[outbound_interface_id].message.payload)

    # Assign the "compiled" update rule as an anomynous function to the schedule entry execute field
    schedule_entry.execute = ( () -> vmp!(node, Val{outbound_interface_id}, rule_arguments...) )

    return schedule_entry
end
