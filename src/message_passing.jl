export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function execute(schedule_entry::ScheduleEntry)
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    isdefined(schedule_entry, :execute) || error("Execute function not defined for schedule entry; perhaps the schedule is not prepared?")
    outbound_dist = schedule_entry.execute()

    # Print output for debugging
    if verbose
        show(schedule_entry)
        println("Result: $(outbound_dist)")
    end

    return outbound_dist
end

# Execute schedules
function execute(schedule::Schedule)
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")

    # Print table header for execution log
    if verbose
        println("\n\nExecution log")
        println("--------------------------------------------")
    end

    for i=1:length(schedule)
        !verbose || println("$(i).")
        execute(schedule[i])
    end

    # Return the last message in the schedule
    entry = schedule[end]
    return entry.node.interfaces[entry.outbound_interface_id].message
end

function clearMessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        clearMessage!(interface)
    end
end

function clearMessages!(edge::Edge)
    # Clear all messages on an edge.
    clearMessage!(edge.head.message)
    clearMessage!(edge.tail.message)
end
