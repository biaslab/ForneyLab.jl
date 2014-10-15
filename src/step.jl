export setReadBuffer, setWriteBuffer, clearBuffers!, addTimeWrap, clearTimeWraps!

setReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=getCurrentGraph()) = (graph.read_buffers[node] = buffer)

setWriteBuffer(interface::Interface, buffer::Vector, graph::FactorGraph=getCurrentGraph()) = (graph.write_buffers[interface] = buffer) # Write buffer for message
setWriteBuffer(edge::Edge, buffer::Vector, graph::FactorGraph=getCurrentGraph()) = (graph.write_buffers[edge] = buffer) # Write buffer for marginal

function clearBuffers!(graph::FactorGraph=getCurrentGraph())
    graph.read_buffers = Dict{TerminalNode, Vector}()
    graph.write_buffers = Dict{Union(Edge,Interface), Vector}()
    return graph
end

function addTimeWrap(from::TerminalNode, to::TerminalNode, graph::FactorGraph=getCurrentGraph())
    !is(from, to) || error("Cannot create time wrap: from and to must be different nodes")
    # Verify that from and to are not already in a time wrap
    for time_wrap in graph.time_wraps
        !(from in time_wrap) || error("Node $(from) is already in another time wrap")
        !(to in time_wrap) || error("Node $(to) is already in another time wrap")
    end

    push!(graph.time_wraps, (from, to))
end

clearTimeWraps!(graph::FactorGraph=getCurrentGraph()) = (graph.time_wraps = Array((TerminalNode, TerminalNode), 0))