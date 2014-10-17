export setReadBuffer, setWriteBuffer, clearBuffers!, addTimeWrap, clearTimeWraps!

function setReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=getCurrentGraph())
    (node in getNodes(graph)) || error("The specified node is not part of the current or specified graph")
    graph.read_buffers[node] = buffer
end

function setWriteBuffer(interface::Interface, buffer::Vector, graph::FactorGraph=getCurrentGraph())
    (interface.node in getNodes(graph)) || error("The specified interface is not part of the current or specified graph")
    graph.write_buffers[interface] = buffer # Write buffer for message
end

function setWriteBuffer(edge::Edge, buffer::Vector, graph::FactorGraph=getCurrentGraph())
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