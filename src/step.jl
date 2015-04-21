export  setReadBuffer,
        setWriteBuffer,
        clearBuffers,
        setTimeWrap,
        clearTimeWraps,
        execute,
        step

function setReadBuffer(node::TerminalNode, buffer::Vector, graph::FactorGraph=current_graph)
    #(node in nodes(graph)) || error("The specified node is not part of the current or specified graph")
    algorithm.read_buffers[node] = buffer
end

function setReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector, graph::FactorGraph=current_graph)
    # Mini-batch assignment for read buffers.
    # buffer is divided over nodes equally.
    n_nodes = length(nodes)
    n_samples_per_node = int(floor(length(buffer)/length(nodes)))
    n_samples_per_node*n_nodes == length(buffer) || error("Buffer length must a multiple of the mini-batch node array length")
    buffmat = reshape(buffer, n_nodes, n_samples_per_node) # samples for one node are present in the rows of buffmat
    for k in 1:n_nodes
        algorithm.read_buffers[nodes[k]] = vec(buffmat[k,:])
    end
    return algorithm.read_buffers[nodes[end]] # Return last node's buffer
end

function setWriteBuffer(interface::Interface, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=current_graph)
    #(interface.node in nodes(graph)) || error("The specified interface is not part of the current or specified graph")
    algorithm.write_buffers[interface] = buffer # Write buffer for message
end

function setWriteBuffer(edge::Edge, buffer::Vector=Array(ProbabilityDistribution,0), graph::FactorGraph=current_graph)
    #(edge in edges(graph)) || error("The specified edge is not part of the current or specified graph")
    algorithm.write_buffers[edge] = buffer # Write buffer for marginal
end

function clearBuffers(graph::FactorGraph=current_graph)
    algorithm.read_buffers = Dict{TerminalNode, Vector}()
    algorithm.write_buffers = Dict{Union(Edge,Interface), Vector}()
    return algorithm
end

function setTimeWrap(from::TerminalNode, to::TerminalNode, storage_graph::FactorGraph=current_graph)
    !is(from, to) || error("Cannot create time wrap: from and to must be different nodes")
    push!(storage_algorithm.time_wraps, (from, to))
end

clearTimeWraps(graph::FactorGraph=current_graph) = (algorithm.time_wraps = Array((TerminalNode, TerminalNode), 0))

function execute(algorithm::Algorithm, graph::FactorGraph=current_graph)
    global current_algorithm = algorithm
    if algorithm.initialized==false
        algorithm.initialize(algorithm.fields)
    end

    return algorithm.execute(algorithm.fields)
end

function step(algorithm::Algorithm, graph::FactorGraph=current_graph)
    # Read buffers
    for (terminal_node, read_buffer) in graph.read_buffers
        !isempty(read_buffer) || error("Read buffer for node $(terminal_node) is empty")
        terminal_node.value = shift!(read_buffer) # pick the first element off the read_buffer
    end

    # Execute schedule
    execute(algorithm, graph)

    # Write buffers
    for (component, write_buffer) in graph.write_buffers
        if typeof(component) == Interface
            isdefined(component.message, :payload) || error("Cannot write message payload to buffer since there is no message present on $(component)")
            push!(write_buffer, deepcopy(component.message.payload))
        elseif typeof(component) == Edge
            push!(write_buffer, calculateMarginal(component))
        end
    end

    # Time wraps
    for (from, to) in graph.time_wraps
        isdefined(from.out.partner.message, :payload) || error("There is no message to move to $(to) for the next timestep")
        to.value = deepcopy(from.out.partner.message.payload)
    end
end