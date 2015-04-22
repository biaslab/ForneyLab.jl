export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function execute(schedule_entry::ScheduleEntry)
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    outbound_interface = schedule_entry.interface
    # Preprocessing: collect all inbound messages and build the inbound_array
    node = outbound_interface.node

    if schedule_entry.message_calculation_rule == sumProduct!
        (outbound_interface_id, inbounds) = SumProduct.collectInbounds(outbound_interface)
    elseif schedule_entry.message_calculation_rule == vmp!
        (outbound_interface_id, inbounds) = VMP.collectInbounds(outbound_interface)
    else
        error("Unknown message calculation rule: $(schedule_entry.message_calculation_rule)")
    end

    # Evaluate message calculation rule
    (rule, outbound_message) = schedule_entry.message_calculation_rule(node, outbound_interface_id, inbounds...)

    # Post processing?
    if isdefined(schedule_entry, :post_processing)
        outbound_message = node.interfaces[outbound_interface_id].message = Message(schedule_entry.post_processing(outbound_message.payload))
    end

    # Print output for debugging
    if verbose && rule != :empty # Internal composite node calls to execute return :empty rule
        interface_name = (name(outbound_interface)!="") ? "$(name(outbound_interface))" : "$(outbound_interface_id)"
        postproc = (isdefined(schedule_entry, :post_processing)) ? string(schedule_entry.post_processing) : ""
        rule_field = "$(rule) $(postproc)"
        println("|$(pad(node.name, 15))|$(pad(interface_name,10))|$(pad(rule_field,30))|$(pad(format(outbound_message.payload),71))|")
    end

    return outbound_message
end

# Execute schedules
function execute(schedule::Schedule)
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")

    # Print table header for execution log
    if verbose
        println("\n|     node      |interface |             rule             |                           calculated message                          |")
        println("|---------------|----------|------------------------------|-----------------------------------------------------------------------|")
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
        interface.message = nothing
    end
end

function clearMessages!(edge::Edge)
    # Clear all messages on an edge.
    edge.head.message = nothing
    edge.tail.message = nothing
end
