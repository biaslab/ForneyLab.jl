function generateScheduleByDFS!(outbound_interface::Interface, 
                                backtrace::Array{Interface, 1}=Array(Interface, 0), 
                                call_list::Array{Interface, 1}=Array(Interface, 0); 
                                allowed_edges=false)
    # Private function to generate a sum product schedule by doing a DFS through the graph.
    # The graph is passed implicitly through the outbound_interface.
    #
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    #
    # outbound_interface: find a schedule to calculate the outbound message on this interface
    # backtrace: backtrace for recursive implementation of DFS
    # call_list: holds the recursive calls
    # allowed_edges: either false or Set{Edge}. If a set is passed, the search will be restricted to edges in this set.
    #
    # Returns: Array{Interface, 1} (not an actual Schedule yet)

    node = outbound_interface.node

    # Apply stopping condition for recursion. When the same interface is called twice, this is indicative of an unbroken loop.
    if outbound_interface in call_list
        # Notify the user to break the loop with an initial message
        error("Loop detected around $(outbound_interface) Consider setting an initial message somewhere in this loop.")
    elseif outbound_interface in backtrace
        # This outbound_interface is already in the schedule
        return backtrace
    else # Stopping condition not satisfied
        push!(call_list, outbound_interface)
    end

    # Check all inbound messages on the other interfaces of the node
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if is(interface, outbound_interface)
            outbound_interface_id = interface_id
        end

        (outbound_interface_id != interface_id) || continue

        if (typeof(allowed_edges)==Set{Edge}) && !(interface.edge in allowed_edges)
            continue
        end

        (interface.partner != nothing) || error("Disconnected interface should be connected: interface #$(interface_id) of $(typeof(node)) $(node.id)")

        if interface.partner.message == nothing # Required message missing.
            if !(interface.partner in backtrace) # Don't recalculate stuff that's already in the schedule.
                # Recursive call
                generateScheduleByDFS!(interface.partner, backtrace, call_list, allowed_edges=allowed_edges)
            end
        end
    end

    # Update call_list and backtrace
    pop!(call_list)

    return push!(backtrace, outbound_interface)
end

function generateSchedule(outbound_interface::Interface; args...)
    # Generate a sum-product Schedule that can be executed to calculate the outbound message on outbound_interface.
    #
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    # The same graph with different messages being present can (and probably will) result in a different schedule.

    return convert(Schedule, generateScheduleByDFS!(outbound_interface; args...))
end

function generateSchedule(partial_schedule::Schedule; args...)
    # Generate a complete schedule based on partial_schedule.
    # A partial schedule only defines the order of a subset of all required messages.
    # This function will find a valid complete schedule that satisfies the partial schedule.
    #
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.

    interface_list = Array(Interface, 0)
    for schedule_entry in partial_schedule
        interface_list = generateScheduleByDFS!(schedule_entry.interface, interface_list; args...)
    end

    return convert(Schedule, interface_list)
end

generateSchedule(partial_list::Array{Interface, 1}; args...) = generateSchedule(convert(Schedule, partial_list); args...)

function generateSchedule(graph::FactorGraph=currentGraph(); args...)
    # Build a sumproduct schedule to calculate all messages towards timewraps and writebuffers
    partial_list = Interface[]
    
    # Collect timewrap interfaces
    for (from_node, to_node) in wraps(graph)
        push!(partial_list, from_node.interfaces[1].partner)
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

    return generateSchedule(partial_list; args...)
end