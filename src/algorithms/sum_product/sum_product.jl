export SumProduct

include("scheduler.jl")

type SumProduct <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
end


############################################
# SumProduct algorithm constructors
############################################

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


############################################
# Type inference and preparation 
############################################

function inferDistributionTypes!(algo::SumProduct)
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}()

    for entry in algo.schedule
        # Generate array of inbound types
        node = entry.node
        outbound_interface_id = entry.outbound_interface_id
        outbound_interface = node.interfaces[outbound_interface_id]

        collectInboundTypes!(entry, node, outbound_interface_id, schedule_entries, algo) # SumProduct specific collection of inbound types
        inferOutboundType!(entry, node, outbound_interface_id, entry.inbound_types, [sumProduct!]) # For the SumProduct algorithm, the only allowed update rule is the sumProduct! rule

        schedule_entries[outbound_interface] = entry # Assign schedule entry to lookup dictionary
    end

    return algo
end

function collectInboundTypes!(entry::ScheduleEntry, node::Node, outbound_interface_id::Int64, schedule_entries::Dict{Interface, ScheduleEntry}, ::SumProduct)
    # Infers inbound types for node relative to the outbound interface
    entry.inbound_types = []

    for (id, interface) in enumerate(node.interfaces)
        if id == outbound_interface_id
            push!(entry.inbound_types, Void) # Outbound interface, push Void
        else
            if interface.partner.message == nothing
                push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
            else # A breaker message is pre-set on the partner interface, push message type
                push!(entry.inbound_types, typeof(interface.partner.message))
            end
        end
    end

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
    for (id, interface) in enumerate(node.interfaces)
        if id == outbound_interface_id
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