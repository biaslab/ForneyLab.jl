export SumProduct

include("scheduler.jl")

type SumProduct <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
end

function SumProduct(graph::FactorGraph=currentGraph())
    # Generates a SumProduct algorithm that propagates messages to all wraps and write buffers.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(graph)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(outbound_interface::Interface)
    # Generates a SumProduct algorithm to calculate the outbound message on outbound_interface.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(outbound_interface)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(partial_list::Vector{Interface})
    # Generates a SumProduct algorithm that at least propagates to all interfaces in the argument vector.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(partial_list)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function SumProduct(edge::Edge)
    # Generates a SumProduct algorithm to calculate the marginal on edge
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule([edge.head, edge.tail])
    function exec(algorithm)
        execute(algorithm.schedule)
        calculateMarginal!(edge)
    end

    algo = SumProduct(exec, schedule)
    inferDistributionTypes!(algo)

    return algo
end

function inferDistributionTypes!(algo::SumProduct)
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}()

    for entry in algo.schedule
        # Generate array of inbound types
        node = entry.node
        outbound_interface_id = entry.outbound_interface_id
        outbound_interface = node.interfaces[outbound_interface_id]

        inferInboundTypes!(entry, node, outbound_interface_id, schedule_entries, algo) # SumProduct specific collection of inbound types
        inferOutboundType!(entry, node, outbound_interface_id, entry.inbound_types, [sumProduct!]) # For the SumProduct algorithm, the only allowed update rule is the sumProduct! rule

        schedule_entries[outbound_interface] = entry # Assign schedule entry to lookup dictionary
    end

    return algo
end

function inferInboundTypes!(entry::ScheduleEntry, node::Node, outbound_interface_id::Int64, schedule_entries::Dict{Interface, ScheduleEntry}, ::SumProduct)
    # Infers inbound types for node relative to the outbound interface
    entry.inbound_types = []
    for i = 1:length(node.interfaces)
        if i == outbound_interface_id
            push!(entry.inbound_types, Void) # Outbouns interface, push Void
        else
            interface = node.interfaces[i]
            if interface.partner.message == nothing
                push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
            else # A breaker message is pre-set on the partner interface, push message type
                push!(entry.inbound_types, typeof(interface.partner.message))
            end
        end
    end

    return entry
end

function inferOutboundType!(entry::ScheduleEntry, node::Node, outbound_interface_id::Int64, inbound_types::Vector, allowed_rules::Vector{Function})
    # Find all compatible calculation rules for the SumProduct algorithm
    outbound_types = []
    for update_function in allowed_rules
        available_rules = methods(update_function, [typeof(node); Type{Val{outbound_interface_id}}; inbound_types; Any])
        for rule in available_rules
            push!(outbound_types, rule.sig.types[end])
        end
    end

    # The outbound outbound_types should contain just one element (there should be just one available)
    if length(outbound_types) == 0
        error("No calculation rule available for inbound types $(inbound_types) on node $(node)")
    elseif length(outbound_types) > 1
        error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types) on node $(node)")
    elseif outbound_types[1] == Any
        # The computation rule produces Any, which indicates that the node is parametrized by its outbound type (e.g. a TerminalNode)
        (typeof(node).parameters[1] <: ProbabilityDistribution) || error("$(typeof(node)) $(node.id) must be parametrized by a ProbabilityDistribution")
        entry.outbound_type = typeof(node).parameters[1]
    elseif outbound_types[1] <: ProbabilityDistribution
        # There is only one possible outbound type and it is a probability distribution
        entry.outbound_type = outbound_types[1]
    else
        error("Unknown output of message calculation rule: $(outbound_types[1]) for node $(node)")
    end

    return entry
end

function inferOutboundType!(entry::ScheduleEntry, node::CompositeNode, outbound_interface_id::Int64, ::Vector, ::Vector{Function})
    # Infer outbound type of composite node
    outbound_interface = node.interfaces[outbound_interface_id]

    # Check if there is already a computation rule defined for this outbound interface
    if !haskey(node.computation_rules, outbound_interface)
        # Try to automatically generate a sum-product algorithm
        clearMessages!(node.internal_graph)
        internal_outbound_interface = node.interfaceid_to_terminalnode[outbound_interface_id].interfaces[1].partner
        composite_algo = SumProduct(internal_outbound_interface)
        node.computation_rules[outbound_interface] = composite_algo
    end

    # Fetch the algorithm that corresponds with outbound interface (either preset or just constructed)
    composite_algo = node.computation_rules[outbound_interface]
    (typeof(composite_algo) <: SumProduct) || error("Only sum product algorithms are currently supported on composite nodes")

    if !isdefined(composite_algo.schedule[end], :outbound_type) # Are the distribution types already inferred?
        inferDistributionTypes!(composite_algo) # Infer types on algorithm coupled to outbound interface of the composite node
    end

    entry.outbound_type = composite_algo.schedule[end].outbound_type # The last entry in the sum-product schedule holds the composite outbound type for this outbound interface

    return entry
end

function prepare!(algo::SumProduct)
    # Populate the graph with vague messages of the correct types
    for entry in algo.schedule
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule (define schedule_entry.execute)
    compile!(algo.schedule, algo)

    return algo
end

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol(sumProduct!)}}, ::InferenceAlgorithm)
    # Generate schedule_entry.execute for schedule entry with sumProduct! update rule

    # Collect references to all required inbound messages for executing message computation rule
    node = schedule_entry.node
    outbound_interface_id = schedule_entry.outbound_interface_id

    rule_arguments = []
    # Add inbound messages to rule_arguments
    for j = 1:length(node.interfaces)
        interface = node.interfaces[j]
        if j == outbound_interface_id
            # Inbound on outbound_interface is irrelevant
            push!(rule_arguments, nothing)
        else
            push!(rule_arguments, interface.partner.message)
        end
    end
    # Add outbound distribution to rule_arguments
    push!(rule_arguments, node.interfaces[outbound_interface_id].message.payload)

    # Assign the "compiled" computation rule as an anomynous function to schedule entry.execute
    schedule_entry.execute = ( () -> sumProduct!(node, Val{outbound_interface_id}, rule_arguments...) )

    return schedule_entry
end