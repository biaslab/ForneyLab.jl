export setReadBuffer, setWriteBuffer, clearBuffers!, addTimeWrap, clearTimeWraps!, step

function setReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=getCurrentGraph())
    (node in getNodes(graph)) || error("The specified node is not part of the current or specified graph")
    graph.read_buffers[node] = buffer
end

function setWriteBuffer(interface::Interface, buffer::Vector=Array(MessagePayload,0), graph::FactorGraph=getCurrentGraph())
    (interface.node in getNodes(graph)) || error("The specified interface is not part of the current or specified graph")
    graph.write_buffers[interface] = buffer # Write buffer for message
end

function setWriteBuffer(edge::Edge, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=getCurrentGraph())
    (edge in getEdges(graph)) || error("The specified edge is not part of the current or specified graph")
    graph.write_buffers[edge] = buffer # Write buffer for marginal
end

function clearBuffers!(graph::FactorGraph=getCurrentGraph())
    graph.read_buffers = Dict{TerminalNode, Vector}()
    graph.write_buffers = Dict{Union(Edge,Interface), Vector}()
    return graph
end

function addTimeWrap(from::TerminalNode, to::TerminalNode, graph::FactorGraph=getCurrentGraph())
    (from in getNodes(graph)) || error("The specified 'from' node is not part of the current or specified graph")
    (to in getNodes(graph)) || error("The specified 'to' node is not part of the current or specified graph")
    !is(from, to) || error("Cannot create time wrap: from and to must be different nodes")
    # Verify that from and to are not already in a time wrap
    for time_wrap in graph.time_wraps
        !(from in time_wrap) || error("Node $(from) is already in another time wrap")
        !(to in time_wrap) || error("Node $(to) is already in another time wrap")
    end

    push!(graph.time_wraps, (from, to))
end

clearTimeWraps!(graph::FactorGraph=getCurrentGraph()) = (graph.time_wraps = Array((TerminalNode, TerminalNode), 0))

function step(graph::FactorGraph=getCurrentGraph(); n_iterations::Int64=1)
    # Reset marginals
    setUninformativeMarginals!(graph)
    # Read buffers
    for (terminal_node, read_buffer) in graph.read_buffers
        !isempty(read_buffer) || error("Read buffer for node $(terminal_node) is empty")
        terminal_node.value = shift!(read_buffer) # pick the first element off the read_buffer
    end
    # Execute schedule
    for iteration = 1:n_iterations
        executeSchedule(graph)
    end
    # Write buffers
    for (component, write_buffer) in graph.write_buffers
        if typeof(component) == Interface
            isdefined(component.message, :payload) || error("Cannot write message payload to buffer since there is no message present on $(component)")
            push!(write_buffer, deepcopy(component.message.payload))
        elseif typeof(component) == Edge
            push!(write_buffer, calculateMarginal(component))
        end
    end
    # Time wrap
    for (from, to) in graph.time_wraps
        isdefined(from.out.partner.message, :payload) || error("There is no message to move to $(to) for the next timestep")
        to.value = deepcopy(from.out.partner.message.payload)
    end
end