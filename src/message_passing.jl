export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function execute(schedule_entry::ScheduleEntry)
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    outbound_interface = schedule_entry.node.interfaces[schedule_entry.outbound_interface_id]

    # Evaluate message calculation rule
    isdefined(schedule_entry, :execute) || error("Execute function not defined for schedule entry; maybe you forgot to prepare the schedule?")
    outbound_dist = schedule_entry.execute()

    # Post processing?
    if isdefined(schedule_entry, :post_processing)
        post_processed_output = schedule_entry.post_processing(outbound_dist)
        if (typeof(post_processed_output) <: ProbabilityDistribution) == false
            # Wrap the output in a DeltaDistribution before packing it in a Message
            post_processed_output = DeltaDistribution(post_processed_output)
        end
        outbound_interface.message = Message(post_processed_output) # Assign new post-processed message to interface
        outbound_dist = outbound_interface.message.payload # Reset outbound to post-processed result
    end

    # Print output for debugging
    if verbose
        node = schedule_entry.node
        interface = node.interfaces[schedule_entry.outbound_interface_id]
        interface_handle = (handle(interface)!="") ? "($(handle(interface)))" : ""
        println(replace("$(schedule_entry.rule) on $(typeof(node)) $(interface.node.id) interface $(schedule_entry.outbound_interface_id) $(interface_handle)", "ForneyLab.", ""))
        if isdefined(schedule_entry, :inbound_types) && isdefined(schedule_entry, :outbound_type)
            println(replace("$(schedule_entry.inbound_types) -> Message{$(schedule_entry.outbound_type)}", "ForneyLab.", ""))
        end
        if isdefined(schedule_entry, :post_processing)
            println(replace("Post processing: $(schedule_entry.post_processing)", "ForneyLab.", ""))
        end
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
        if verbose
            println("$(i).")
        end
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