function generateSumProductSchedule(outbound_interface::Interface; args...)
    # Generate a sum-product Schedule that can be executed to calculate the outbound message on outbound_interface.

    return convert(Schedule, generateScheduleByDFS!(outbound_interface; args...), sumProductRule!)
end

function generateSumProductSchedule(partial_schedule::Schedule; args...)
    # Generate a complete schedule based on partial_schedule.
    # A partial schedule only defines the order of a subset of all required messages.
    # This function will find a valid complete schedule that satisfies the partial schedule.

    interface_list = Array(Interface, 0)
    for schedule_entry in partial_schedule
        outbound_interface = schedule_entry.node.interfaces[schedule_entry.outbound_interface_id]
        interface_list = generateScheduleByDFS!(outbound_interface, interface_list; args...)
    end

    return convert(Schedule, interface_list, sumProductRule!)
end

generateSumProductSchedule(partial_list::Array{Interface, 1}; args...) = generateSumProductSchedule(convert(Schedule, partial_list, sumProductRule!); args...)

function generateSumProductSchedule(graph::FactorGraph=currentGraph(); args...)
    # Build a sumproduct schedule to calculate all messages towards wraps and writebuffers
    partial_list = Interface[]

    # Collect wrap interfaces
    for wrap in wraps(graph)
        push!(partial_list, wrap.tail.interfaces[1].partner)
        #push!(partial_list, wrap.head.interfaces[1].partner)
    end

    # Collect write buffer interfaces
    for entry in keys(graph.write_buffers)
        if typeof(entry) == Interface
            push!(partial_list, entry)
        elseif typeof(entry) == Edge
            push!(partial_list, entry.head)
            push!(partial_list, entry.tail)
        end
    end

    return generateSumProductSchedule(partial_list; args...)
end
