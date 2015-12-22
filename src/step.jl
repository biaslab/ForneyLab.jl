export  attachReadBuffer,
        attachWriteBuffer,
        detachReadBuffer,
        detachWriteBuffer,
        detachBuffers,
        emptyWriteBuffers,
        execute,
        step,
        run

import Base.run
import Base.step

function attachReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=currentGraph())
    hasNode(graph, node) || error("The specified node is not part of the current or specified graph")
    graph.read_buffers[node] = buffer
end

function attachReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector, graph::FactorGraph=currentGraph())
    # Mini-batch assignment for read buffers.
    # buffer is divided over nodes equally.
    n_nodes = length(nodes)
    n_samples_per_node = round(Int, floor(length(buffer)/length(nodes)))
    n_samples_per_node*n_nodes == length(buffer) || error("Buffer length must a multiple of the mini-batch node array length")
    buffmat = reshape(buffer, n_nodes, n_samples_per_node) # samples for one node are present in the rows of buffmat
    for k in 1:n_nodes
        hasNode(graph, nodes[k]) || error("One of the specified nodes is not part of the current or specified graph")
        graph.read_buffers[nodes[k]] = vec(buffmat[k,:])
    end

    return graph.read_buffers[nodes[end]] # Return last node's buffer
end

function detachReadBuffer(nd::TerminalNode, graph::FactorGraph=currentGraph())
    hasNode(graph, nd) || error("The specified node is not part of the current or specified graph")
    haskey(graph.read_buffers, nd) || error("There is no read buffer attached to the specified node")

    delete!(graph.read_buffers, nd)
    return graph
end

function attachWriteBuffer(interface::Interface, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=currentGraph())
    hasNode(graph, interface.node) || error("The specified interface is not part of the current or specified graph")
    graph.write_buffers[interface] = buffer # Write buffer for message
end

# function attachWriteBuffer(interfaces::Vector{Interface}, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=currentGraph())
#     # Mini-batch assignment for write buffers.
#     # After each step the batch results are appended to the buffer

# end

function detachWriteBuffer(interface::Interface, graph::FactorGraph=currentGraph())
    hasNode(graph, interface.node) || error("The specified interface is not part of the current or specified graph")
    haskey(graph.write_buffers, interface) || error("There is no write buffer attached to the specified interface")

    delete!(graph.write_buffers, interface)
    return graph
end

function attachWriteBuffer(edge::Edge, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=currentGraph())
    hasEdge(graph, edge) || error("The specified edge is not part of the current or specified graph")
    graph.write_buffers[edge] = buffer # Write buffer for marginal
end

function detachWriteBuffer(edge::Edge, graph::FactorGraph=currentGraph())
    hasEdge(graph, edge) || error("The specified edge is not part of the current or specified graph")
    haskey(graph.write_buffers, edge) || error("There is no write buffer attached to the specified edge")

    delete!(graph.write_buffers, edge)
    return graph
end

function detachBuffers(graph::FactorGraph=currentGraph())
    graph.read_buffers = Dict{TerminalNode, Vector}()
    graph.write_buffers = Dict{Union{Edge,Interface}, Vector}()
end

function emptyWriteBuffers(graph::FactorGraph=currentGraph())
    for (k, v) in graph.write_buffers
        empty!(v) # Clear the vector but keep the pointer
    end
end

function execute(algorithm::InferenceAlgorithm)
    # Call algorithm's execute function with itself as argument
    # prepare!(algorithm) should always be called before the first call to execute(algorithm)

    return algorithm.execute(algorithm)
end

function step(algorithm::InferenceAlgorithm)
    # Execute algorithm for 1 timestep.
    # prepare!(algorithm) should always be called before the first call to step(algorithm)

    # Read buffers
    for (terminal_node, read_buffer) in currentGraph().read_buffers
        !isempty(read_buffer) || error("Read buffer for node $(terminal_node) is empty")
        terminal_node.value = shift!(read_buffer) # pick the first element off the read_buffer
    end

    # Execute schedule
    result = execute(algorithm)

    # Write buffers
    for (component, write_buffer) in currentGraph().write_buffers
        if typeof(component) == Interface
            push!(write_buffer, deepcopy(component.message.payload))
        elseif typeof(component) == Edge
            push!(write_buffer, calculateMarginal(component))
        end
    end

    # Wraps
    for wrap in wraps(currentGraph())
        wrap.sink.value = deepcopy(wrap.source.interfaces[1].partner.message.payload)
    end

    return result
end

function run(algorithm::InferenceAlgorithm)
    # Call step(algorithm) repeatedly until at least one read buffer is exhausted
    prepare!(algorithm)
    if length(currentGraph().read_buffers) > 0
        while !any(isempty, values(currentGraph().read_buffers))
            step(algorithm)
        end
    else
        # No read buffers, just call step once
        step(algorithm)
    end
end
