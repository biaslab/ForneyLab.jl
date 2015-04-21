export generateSchedule!

function generateSchedule!(outbound_interface::Interface, scheme::InferenceScheme=currentScheme(); args...)
    schedule = generateSchedule(outbound_interface, scheme; args...)
    return subgraph(scheme, outbound_interface.edge).internal_schedule = schedule
end

function generateSchedule!(partial_schedule::Schedule, scheme::InferenceScheme=currentScheme(); args...)
    schedule = generateSchedule(partial_schedule, scheme; args...)
    return subgraph(scheme, partial_schedule[1].edge).internal_schedule = schedule
end
generateSchedule!(partial_list::Array{Interface, 1}, scheme::InferenceScheme=currentScheme(); args...) = generateSchedule!(convert(Schedule, partial_list), scheme; args...)

function generateSchedule!(sg::Subgraph, scheme::InferenceScheme=currentScheme())
    # Generate an internal and external schedule for the subgraph

    interface_list_for_univariate = Array(Interface, 0)
    internal_interface_list = Array(Interface, 0)
    sg.internal_schedule = Array(ScheduleEntry, 0)
    # The internal schedule makes sure that incoming internal messages over internal edges connected to nodes (g) are present
    for g_node in nodesConnectedToExternalEdges(sg) # All nodes that are connected to at least one external edge
        outbound_interfaces = Array(Interface, 0) # Array that holds required outbound for the case of one internal edge connected to g_node
        for interface in g_node.interfaces
            if interface.edge in sg.internal_edges # edge carries incoming internal message
                # Store outbound interfaces for check later on
                if !(interface in internal_interface_list) && !(interface in interface_list_for_univariate)
                    push!(outbound_interfaces, interface) # If we were to add the outbound to the schedule (for the case of univariate q), this is the one
                end

                # Extend internal_schedule to calculate the inbound message on interface
                try
                    internal_interface_list = generateScheduleByDFS(interface.partner, internal_interface_list, Array(Interface, 0), scheme, stay_in_subgraph=true)
                catch
                    error("Cannot generate internal schedule for possibly loopy subgraph with internal edge $(interface.edge).")
                end
            end
        end

        # For the case that g_node is connected to one internal edge,
        # the calculation reduces to the naive vmp update which requires the outbound (Dauwels, 2007)
        if length(outbound_interfaces) == 1
            push!(interface_list_for_univariate, outbound_interfaces[1])
        end
    end

    # Make sure that messages are propagated to the timewraps
    interface_list_for_time_wraps = Array(Interface, 0)
    for (from_node, to_node) in scheme.time_wraps
        if subgraph(scheme, from_node.out.edge) == sg # Timewrap is the responsibility of this subgraph
            interface_list_for_time_wraps = [interface_list_for_time_wraps, generateScheduleByDFS(from_node.out.partner, Array(Interface, 0), Array(Interface, 0), scheme, stay_in_subgraph=true)]
        end
    end

    # Make sure that messages are propagated to the write buffers
    interface_list_for_write_buffers = Array(Interface, 0)
    for entry in keys(scheme.write_buffers)
        if typeof(entry) == Interface
            interface_list_for_write_buffers = [interface_list_for_write_buffers, generateScheduleByDFS(entry, Array(Interface, 0), Array(Interface, 0), scheme, stay_in_subgraph=true)]
        elseif typeof(entry) == Edge
            interface_list_for_write_buffers = [interface_list_for_write_buffers, generateScheduleByDFS(entry.head, Array(Interface, 0), Array(Interface, 0), scheme, stay_in_subgraph=true)]
            interface_list_for_write_buffers = [interface_list_for_write_buffers, generateScheduleByDFS(entry.tail, Array(Interface, 0), Array(Interface, 0), scheme, stay_in_subgraph=true)]
        end
    end

    # Schedule for univariate comes after internal schedule, because it can depend on inbounds
    sg.internal_schedule = convert(Schedule, unique([internal_interface_list, interface_list_for_univariate, interface_list_for_time_wraps, interface_list_for_write_buffers]))

    return sg
end

function generateSchedule!(scheme::InferenceScheme=currentScheme())
    for subgraph in scheme.factorization
        generateSchedule!(subgraph, scheme)
    end
    return scheme
end
