export  attachReadBuffer,
        attachWriteBuffer,
        detachReadBuffer,
        detachWriteBuffer,
        detachBuffers,
        emptyWriteBuffers,
        execute,
        step,
        run

function attachReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=current_graph)
    hasNode(graph, node) || error("The specified node is not part of the current or specified graph")
    graph.read_buffers[node] = buffer
end

function attachReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector, graph::FactorGraph=current_graph)
    # Mini-batch assignment for read buffers.
    # buffer is divided over nodes equally.
    # TODO: this function should be renamed to reflect what it does
    n_nodes = length(nodes)
    n_samples_per_node = int(floor(length(buffer)/length(nodes)))
    n_samples_per_node*n_nodes == length(buffer) || error("Buffer length must a multiple of the mini-batch node array length")
    buffmat = reshape(buffer, n_nodes, n_samples_per_node) # samples for one node are present in the rows of buffmat
    for k in 1:n_nodes
        hasNode(graph, nodes[k]) || error("One of the specified nodes is not part of the current or specified graph")
        graph.read_buffers[nodes[k]] = vec(buffmat[k,:])
    end

    return graph.read_buffers[nodes[end]] # Return last node's buffer
end

function detachReadBuffer(nd::TerminalNode, graph::FactorGraph=current_graph)
    hasNode(graph, nd) || error("The specified node is not part of the current or specified graph")
    haskey(graph.read_buffers, nd) || error("There is no read buffer attached to the specified node")

    delete!(graph.read_buffers, nd)
    return graph
end

function attachWriteBuffer(interface::Interface, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=current_graph)
    hasNode(graph, interface.node) || error("The specified interface is not part of the current or specified graph")
    graph.write_buffers[interface] = buffer # Write buffer for message
end

function detachWriteBuffer(interface::Interface, graph::FactorGraph=current_graph)
    hasNode(graph, interface.node) || error("The specified interface is not part of the current or specified graph")
    haskey(graph.write_buffers, interface) || error("There is no write buffer attached to the specified interface")

    delete!(graph.write_buffers, interface)
    return graph
end

function attachWriteBuffer(edge::Edge, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=current_graph)
    hasEdge(graph, edge) || error("The specified edge is not part of the current or specified graph")
    graph.write_buffers[edge] = buffer # Write buffer for marginal
end

function detachWriteBuffer(edge::Edge, graph::FactorGraph=current_graph)
    hasEdge(graph, edge) || error("The specified edge is not part of the current or specified graph")
    haskey(graph.write_buffers, edge) || error("There is no write buffer attached to the specified edge")

    delete!(graph.write_buffers, edge)
    return graph
end

function detachBuffers(graph::FactorGraph=current_graph)
    graph.read_buffers = Dict{TerminalNode, Vector}()
    graph.write_buffers = Dict{Union(Edge,Interface), Vector}()
end

function emptyWriteBuffers(graph::FactorGraph=current_graph)
    for (k, v) in graph.write_buffers
        empty!(v) # Clear the vector but keep the pointer
    end
end

function execute(algorithm::Algorithm, graph::FactorGraph=current_graph)
    # Execute algorithm on graph
    global current_algorithm = algorithm
    return algorithm.execute(algorithm.fields)
end

function step(algorithm::Algorithm, graph::FactorGraph=current_graph)
    # Execute algorithm for 1 timestep. 

    # Read buffers
    for (terminal_node, read_buffer) in graph.read_buffers
        !isempty(read_buffer) || error("Read buffer for node $(terminal_node) is empty")
        terminal_node.value = shift!(read_buffer) # pick the first element off the read_buffer
    end

    # Execute schedule
    result = execute(algorithm, graph)

    # Write buffers
    for (component, write_buffer) in graph.write_buffers
        if typeof(component) == Interface
            isdefined(component.message, :payload) || error("Cannot write message payload to buffer since there is no message present on $(component)")
            push!(write_buffer, deepcopy(component.message.payload))
        elseif typeof(component) == Edge
            push!(write_buffer, calculateMarginal(component))
        end
    end

    # Wraps
    for wrap in wraps(graph)
        isdefined(wrap.source.interfaces[1].partner.message, :payload) || error("There is no message to move to $(wrap.sink) for the next timestep")
        wrap.sink.value = deepcopy(wrap.source.interfaces[1].partner.message.payload)
    end

    return result
end

function run(algorithm::Algorithm, graph::FactorGraph=current_graph)
    # Call step(algorithm, graph) repeatedly until at least one read buffer is exhausted
    if length(graph.read_buffers) > 0
        while !any(isempty, values(graph.read_buffers))
            step(algorithm, graph)
        end
    else
        # No read buffers, just call step once
        step(algorithm, graph)
    end
end