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
    compile!(algo.schedule, algo)

    return algo
end

function SumProduct(outbound_interface::Interface)
    # Generates a SumProduct algorithm to calculate the outbound message on outbound_interface.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(outbound_interface)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(exec, schedule)
    compile!(algo.schedule, algo)

    return algo
end

function SumProduct(partial_list::Vector{Interface})
    # Generates a SumProduct algorithm that at least propagates to all interfaces in the argument vector.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(partial_list)
    exec(algorithm) = execute(algorithm.schedule)

    algo = SumProduct(exec, schedule)
    compile!(algo.schedule, algo)

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
    compile!(algo.schedule, algo)

    return algo
end

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol("ForneyLab.sumProduct!")}}, ::InferenceAlgorithm)
    # Compile ScheduleEntry objects for SumProduct algorithm
    # Generates schedule_entry.execute function

    # Collect references to all required inbound messages for executing message computation rule
    outbound_interface = schedule_entry.interface
    node = outbound_interface.node
    outbound_interface_index = 0
    inbounds = Any[]
    for j = 1:length(node.interfaces)
        interface = node.interfaces[j]
        if is(interface, outbound_interface)
            outbound_interface_index = j
            push!(inbounds, nothing)
            continue
        end
        push!(inbounds, interface.partner.message)
    end

    # Assign the "compiled" update rule as an anomynous function to the schedule entry execute field
    schedule_entry.execute = ( () -> sumProduct!(node, outbound_interface_index, inbounds...) )

    return schedule_entry
end
