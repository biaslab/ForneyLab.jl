export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function execute(schedule_entry::ScheduleEntry)
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    outbound_interface = schedule_entry.interface
    node = outbound_interface.node

    # Evaluate message calculation rule
    (rule, outbound_message) = schedule_entry.execute()

    # Post processing?
    if isdefined(schedule_entry, :post_processing)
        post_processed_output = schedule_entry.post_processing(outbound_message.payload)
        if (typeof(post_processed_output) <: ProbabilityDistribution) == false
            # Wrap the output in a DeltaDistribution before packing it in a Message
            post_processed_output = DeltaDistribution(post_processed_output)
        end
        outbound_message = node.interfaces[outbound_interface_index].message = Message(post_processed_output)
    end

    # Print output for debugging
    if verbose && rule != :empty
        interface_handle = (handle(outbound_interface)!="") ? "$(handle(outbound_interface))" : "$(outbound_interface_index)"
        postproc = (isdefined(schedule_entry, :post_processing)) ? string(schedule_entry.post_processing) : ""
        rule_field = "$(rule) $(postproc)"
        println("$(node.id) [$(interface_handle)], $(rule_field): $(format(outbound_message.payload))")
    end

    return outbound_message
end

# Execute schedules
function execute(schedule::Schedule)
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")

    # Print table header for execution log
    if verbose
        println("\n\nExecution log (node [interface], rule: result)")
        println("--------------------------------------------")
    end

    for schedule_entry in schedule
        execute(schedule_entry)
    end
    # Return the last message in the schedule
    return schedule[end].interface.message
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
