function generateSumProductSchedule(outbound_interface::Interface; args...)
    # Generate a sum-product Schedule that can be executed to calculate the outbound message on outbound_interface.

    return convert(Schedule, generateScheduleByDFS!(outbound_interface; args...), SumProductRule)
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

    return convert(Schedule, interface_list, SumProductRule)
end

generateSumProductSchedule(partial_list::Array{Interface, 1}; args...) = generateSumProductSchedule(convert(Schedule, partial_list, SumProductRule); args...)

generateSumProductSchedule(graph::FactorGraph=currentGraph(); args...) = generateSumProductSchedule(interfacesFacingWrapsOrBuffers(graph); args...)
