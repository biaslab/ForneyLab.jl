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
        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        inbound_types = []
        for i=1:length(entry.node.interfaces)
            if i == entry.outbound_interface_id
                push!(inbound_types, Void)
            else
                interface = entry.node.interfaces[i]
                push!(inbound_types, schedule_entries[interface].outbound_type)
            end
        end

        # Find all compatible calculation rules
        available_rules = methods(sumProduct!, [typeof(entry.node), Type{Val{entry.outbound_interface_id}}, inbound_types, Any])
        outbound_types = [rule.sig.types[end].parameters[1] for rule in available_rules]

        # The oubound outbound_types should contain just one element (there should be just one available)
        if length(outbound_types) == 0
            error("No calculation rule available for inbound types $(inbound_types).\n$(entry)")
        elseif length(outbound_types) > 1
            error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types).\n$(entry)")
        elseif outbound_types[1] == Any
            # The computation rule produces Any, which indicates that the node is parametrized by its outbound type (e.g. a TerminalNode)
            (entry.node.parameters[1] <: ProbabilityDistribution) || error("$(typeof(entry.node)) $(entry.node.id) must be parametrized by a ProbabilityDistribution")
            entry.outbound_type = entry.node.parameters[1]
        else
            entry.outbound_type = outbound_types[1]
        end

        schedule_entries[outbound_interface] = entry # Assign schedule entry to lookup dictionary
    end

    return algo
end

function prepare!(algo::SumProduct)
    # Populate the graph with vague messages of the correct types
    for entry in algo.schedule
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule (define schedule_entry.execute)
    compile!(algo.schedule)

    return algo
end

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol("ForneyLab.sumProduct!")}}, ::InferenceAlgorithm)
    # Generate schedule_entry.execute for SumProduct algorithm

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
    push!(rule_arguments, interface.message.payload)

    # Assign the "compiled" computation rule as an anomynous function to schedule entry.execute
    schedule_entry.execute = ( () -> sumProduct!(node, Val{outbound_interface_id}, rule_arguments...) )

    return schedule_entry
end
