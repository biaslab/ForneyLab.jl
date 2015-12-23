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

function VariationalBayes(graph::FactorGraph=currentGraph(); n_iterations::Int64=50)
    # Generates a VariationalBayes algorithm that propagates messages to all write buffers and wraps.
    # Uses a mean field factorization and autoscheduler
    factorization = factorize(graph) # Mean field factorization
    generateVariationalBayesSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = vagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        vagueQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
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
    q_distributions = vagueQDistributions(factorization) # Initialize vague q distributions

    function exec(algorithm)
        vagueQDistributions!(algorithm.q_distributions) # Reset q distributions before next step
        for iteration = 1:algorithm.n_iterations
            execute(algorithm.factorization, algorithm.q_distributions) # For all subgraphs, execute internal and external schedules
        end
    end

    algo = VariationalBayes(exec, factorization, q_distributions, n_iterations)
    inferDistributionTypes!(algo)

    return algo
end
VariationalBayes(cluster_edges...; n_iterations::Int64=50) = VariationalBayes(currentGraph(), cluster_edges...; n_iterations=n_iterations)

function inferDistributionTypes!(algo::VariationalBayes)
    # Infer the payload types for all messages in the internal schedules

    for factor in algo.factorization.factors
        # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
        schedule = factor.internal_schedule
        schedule_entries = Dict{Interface, ScheduleEntry}()

        for entry in schedule
            # Generate array of inbound types
            outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
            inbound_types = []

            # Collect references to all required inbound messages for executing message computation rule
            for i = 1:length(entry.node.interfaces)
                if i == entry.outbound_interface_id
                    push!(inbound_types, Void) # This interface is outbound
                else
                    interface = entry.node.interfaces[i]
                    if !haskey(algo.factorization.edge_to_subgraph, interface.edge) || !haskey(algo.factorization.edge_to_subgraph, outbound_interface.edge)
                        # Inbound and/or outbound edge is not explicitly listed in the factorization edge list.
                        # This is possible if one of those edges is internal to a composite node.
                        # We will default to sum-product message passing, and consume the message on the inbound interface.
                        # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph,
                        # so it is safe to assume sumproduct.
                        try push!(inbound_types, Message{schedule_entries[interface.partner].outbound_type}) catch error("$(interface) is not connected to an edge.") end
                        break
                    end

                    # Should we require the inbound message or marginal?
                    if is(algo.factorization.edge_to_subgraph[interface.edge], algo.factorization.edge_to_subgraph[outbound_interface.edge])
                        # Both edges in same subgraph, require message
                        try push!(inbound_types, Message{schedule_entries[interface.partner].outbound_type}) catch error("$(interface) is not connected to an edge.") end
                    else
                        # A subgraph border is crossed, require marginal
                        # The factor is the set of internal edges that are in the same subgraph
                        sg = algo.factorization.edge_to_subgraph[interface.edge]
                        push!(inbound_types, typeof(algo.q_distributions[(outbound_interface.node, sg)].distribution))
                    end
                end
            end

            # Find all compatible calculation rules
            available_sumproduct_rules = methods(sumProduct!, [typeof(entry.node); Type{Val{entry.outbound_interface_id}}; inbound_types; Any])
            outbound_sumproduct_types = [rule.sig.types[end] for rule in available_sumproduct_rules]

            available_vmp_rules = methods(vmp!, [typeof(entry.node); Type{Val{entry.outbound_interface_id}}; inbound_types; Any])
            outbound_vmp_types = [rule.sig.types[end] for rule in available_vmp_rules]

            outbound_types = [outbound_sumproduct_types; outbound_vmp_types]

            # Assign the outbound type to the schedule entry
            # The outbound outbound_types should contain just one element (there should be just one available)
            if length(outbound_types) == 0
                error("No calculation rule available for inbound types $(inbound_types).\n$(entry)")
            elseif length(outbound_types) > 1
                error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types).\n$(entry)")
            elseif outbound_types[1] == Any
                # The computation rule produces Any, which indicates that the node is parametrized by its outbound type (e.g. a TerminalNode)
                (typeof(entry.node).parameters[1] <: ProbabilityDistribution) || error("$(typeof(entry.node)) $(entry.node.id) must be parametrized by a ProbabilityDistribution")
                entry.outbound_type = typeof(entry.node).parameters[1]
            elseif outbound_types[1] <: ProbabilityDistribution
                entry.outbound_type = outbound_types[1]
            else
                error("Unknown output of message calculation rule: $(outbound_types[1])\n$(entry)")
            end
        end
    end

    return algo
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
    # Generate schedule_entry.execute for vmp update

    # Collect references to all required inbound messages for executing message computation rule
    node = schedule_entry.node
    outbound_interface_id = schedule_entry.outbound_interface_id

    rule_arguments = []
    # Add inbound messages to rule_arguments
    for j = 1:length(node.interfaces)
        interface = node.interfaces[j]
        if j == outbound_interface_id
            # Inbound on outbound_interface is irrelevant
            push!(rule_arguments, nothing) # This interface is outbound, push "nothing"
        else
            if !haskey(algo.factorization.edge_to_subgraph, interface.edge) || !haskey(algo.factorization.edge_to_subgraph, outbound_interface.edge)
                # Inbound and/or outbound edge is not explicitly listed in the algo.fields.
                # This is possible if one of those edges is internal to a composite node.
                # We will default to sum-product message passing, and consume the message on the inbound interface.
                # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
                try push!(rule_arguments, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
                break
            end

            # Should we require the inbound message or marginal?
            if is(algo.factorization.edge_to_subgraph[interface.edge], algo.factorization.edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                try push!(rule_arguments, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = algo.factorization.edge_to_subgraph[interface.edge]
                push!(rule_arguments, algo.q_distributions[(outbound_interface.node, sg)].distribution)
            end
        end
    end
    # Add outbound distribution to rule_arguments
    push!(rule_arguments, node.interfaces[outbound_interface_id].message.payload)

    # Assign the "compiled" update rule as an anomynous function to the schedule entry execute field
    schedule_entry.execute = ( () -> vmp!(node, outbound_interface_index, rule_arguments...) )

    return schedule_entry
end
