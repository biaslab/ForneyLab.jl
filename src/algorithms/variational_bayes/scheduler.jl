# Generates schedules for the VMP algorithm.
# VMP schedules (internal and external) are stored on the subgraphs in the factorization.
# There are no call signatures for message passing to specific interfaces or edges (as with sum product);
# when required, these should be set through write buffers or wraps.

function generateVariationalBayesSchedule!(f::RecognitionFactorization, graph::FactorGraph=currentGraph())
    # Generate and store an internal and external schedule for each subgraph in the factorization
    for sg in f.factors
        generateVariationalBayesSchedule!(sg, graph)
    end
    return f
end

function generateVariationalBayesSchedule!(sg::Subgraph, graph::FactorGraph=currentGraph())
    # Generate and store an internal and external schedule for the subgraph

    external_schedule = nodesConnectedToExternalEdges(sg) # All nodes that are connected to at least one external edge
    interface_list_for_univariate = Array(Interface, 0) # This array holds the interfaces required to calculate univariate marginals (Dauwels, 2007)
    internal_interface_list = Array(Interface, 0) # This array will hold the sumproduct schedule for internal message passing
    external_nodes_interfaces = Array(Interface, 0) # This array will hold the interfaces that need to be calculated with a variationalRule! update
    sg.internal_schedule = Array(ScheduleEntry, 0)
    # The internal schedule makes sure that incoming internal messages over internal edges connected to nodes (g) are present
    for g_node in external_schedule # Internal schedule depends on nodes connected to external edges
        outbound_interfaces = Array(Interface, 0) # Array that holds required outbound for the case of one internal edge connected to g_node
        for interface in g_node.interfaces
            if interface.edge in sg.internal_edges # edge carries incoming internal message
                # Store outbound interfaces in care we need to calculate an univariate marginal later on
                if !(interface in internal_interface_list) && !(interface in interface_list_for_univariate)
                    push!(outbound_interfaces, interface) # If we were to eventually add the outbound to the schedule (for the case of univariate q), this is the one if it is not already in the interface list
                end

                try # Extend internal_interface_list to calculate the inbound message on interface
                    internal_interface_list = generateScheduleByDFS!(interface.partner, internal_interface_list, Array(Interface, 0), allowed_edges=sg.internal_edges)
                catch
                    error("Cannot generate internal schedule for possibly loopy subgraph with internal edge $(interface.edge).")
                end
            end
        end

        # For the case that g_node is connected to one internal edge,
        # the calculation reduces to the naive vmp update (Dauwels, 2007). This update requires the outbound that was stored earlier.
        if length(outbound_interfaces) == 1
            push!(interface_list_for_univariate, outbound_interfaces[1])
        end

        # Store the interfaces connected to external nodes. These need to be calculated with a vmp update rule instead of sumproduct.
        external_nodes_interfaces = [external_nodes_interfaces; g_node.interfaces]
    end

    # Make sure that messages are propagated to the wraps
    interface_list_for_wraps = Array(Interface, 0)
    sg_nodes = nodes(sg)
    for wrap in wraps(graph)
        if wrap.source in sg_nodes # Timewrap is the responsibility of this subgraph
            interface_list_for_wraps = [interface_list_for_wraps; generateScheduleByDFS!(wrap.source.interfaces[1].partner, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
        end
    end

    # Make sure that messages are propagated to the write buffers
    interface_list_for_write_buffers = Array(Interface, 0)
    for entry in keys(graph.write_buffers)
        if typeof(entry) == Interface
            if entry.edge in sg.internal_edges # Interface is the responsibility of sg
                interface_list_for_write_buffers = [interface_list_for_write_buffers; generateScheduleByDFS!(entry, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
            end
        elseif typeof(entry) == Edge
            if entry in sg.internal_edges # Edge is the responsibility of sg
                interface_list_for_write_buffers = [interface_list_for_write_buffers; generateScheduleByDFS!(entry.head, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
                interface_list_for_write_buffers = [interface_list_for_write_buffers; generateScheduleByDFS!(entry.tail, Array(Interface, 0), Array(Interface, 0), allowed_edges=sg.internal_edges)]
            end
        end
    end

    # Convert the interface list to a schedule. The schedule for univariate comes after the internal schedule, because it can depend on inbound messages calculated earlier
    schedule = convert(Schedule, unique([internal_interface_list; interface_list_for_univariate; interface_list_for_wraps; interface_list_for_write_buffers]), sumProductRule!)
    # Convert to a vmp update when a schedule entry interface belongs to an external node
    for entry in schedule
        if entry.node.interfaces[entry.outbound_interface_id] in external_nodes_interfaces
            entry.rule = variationalRule!
        end
    end

    # Store schedules
    sg.internal_schedule = schedule
    sg.external_schedule = collect(external_schedule) # Convert set to array

    return sg
end
