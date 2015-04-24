# Generates schedules for the VMP algorithm.
# VMP schedules (internal and external) are stored on the subgraphs in the factorization.
# There are no call signatures for message passing to specific interfaces or edges;
# when required, these should be set through write buffers or time wraps.

function generateSchedule!(f::Factorization, graph::FactorGraph=current_graph)
    # Generate and store an internal and external schedule for each subgraph in the factorization
    for sg in f.factors
        generateSchedule!(sg, graph)
    end
    return f
end

function generateSchedule!(sg::Subgraph, graph::FactorGraph=current_graph)
    # Generate and store an internal and external schedule for the subgraph

    external_schedule = nodesConnectedToExternalEdges(sg) # All nodes that are connected to at least one external edge

    interface_list_for_univariate = Array(Interface, 0)
    internal_interface_list = Array(Interface, 0)
    sg.internal_schedule = Array(ScheduleEntry, 0)
    # The internal schedule makes sure that incoming internal messages over internal edges connected to nodes (g) are present
    for g_node in external_schedule # Internal schedule depends on nodes connected to external edges
        outbound_interfaces = Array(Interface, 0) # Array that holds required outbound for the case of one internal edge connected to g_node
        for interface in g_node.interfaces
            if interface.edge in sg.internal_edges # edge carries incoming internal message
                # Store outbound interfaces for check later on
                if !(interface in internal_interface_list) && !(interface in interface_list_for_univariate)
                    push!(outbound_interfaces, interface) # If we were to add the outbound to the schedule (for the case of univariate q), this is the one
                end

                # Extend internal_schedule to calculate the inbound message on interface
                try
                    internal_interface_list = SumProduct.generateScheduleByDFS(interface.partner, internal_interface_list, Array(Interface, 0), allowed_edges=sg.internal_edges)
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
    for (from_node, to_node) in graph.time_wraps
        if from_node.out.edge in sg.internal_edges # Timewrap is the responsibility of this subgraph
            interface_list_for_time_wraps = [interface_list_for_time_wraps, SumProduct.generateScheduleByDFS(from_node.out.partner, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
        end
    end

    # Make sure that messages are propagated to the write buffers
    interface_list_for_write_buffers = Array(Interface, 0)
    for entry in keys(graph.write_buffers)
        if typeof(entry) == Interface
            interface_list_for_write_buffers = [interface_list_for_write_buffers, SumProduct.generateScheduleByDFS(entry, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
        elseif typeof(entry) == Edge
            interface_list_for_write_buffers = [interface_list_for_write_buffers, SumProduct.generateScheduleByDFS(entry.head, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
            interface_list_for_write_buffers = [interface_list_for_write_buffers, SumProduct.generateScheduleByDFS(entry.tail, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
        end
    end

    # Schedule for univariate comes after internal schedule, because it can depend on inbounds
    sg.internal_schedule = convert(Schedule, unique([internal_interface_list, interface_list_for_univariate, interface_list_for_time_wraps, interface_list_for_write_buffers]))
    sg.external_schedule = external_schedule

    return sg 
end