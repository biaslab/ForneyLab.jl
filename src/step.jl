export  setReadBuffer,
        setWriteBuffer,
        clearBuffers,
        setTimeWrap,
        clearTimeWraps,
        step

function setReadBuffer(node::TerminalNode, buffer::Vector, scheme::InferenceScheme=currentGraph().active_scheme)
    #(node in nodes(graph)) || error("The specified node is not part of the current or specified graph")
    scheme.read_buffers[node] = buffer
end

function setReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector, scheme::InferenceScheme=currentGraph().active_scheme)
    # Mini-batch assignment for read buffers.
    # buffer is divided over nodes equally.
    n_nodes = length(nodes)
    n_samples_per_node = int(floor(length(buffer)/length(nodes)))
    n_samples_per_node*n_nodes == length(buffer) || error("Buffer length must a multiple of the mini-batch node array length")
    buffmat = reshape(buffer, n_nodes, n_samples_per_node) # samples for one node are present in the rows of buffmat
    for k in 1:n_nodes
        scheme.read_buffers[nodes[k]] = vec(buffmat[k,:])
    end
    return scheme.read_buffers[nodes[end]] # Return last node's buffer
end

function setWriteBuffer(interface::Interface, buffer::Vector=Array(ProbabilityDistribution,0), scheme::InferenceScheme=currentGraph().active_scheme)
    #(interface.node in nodes(graph)) || error("The specified interface is not part of the current or specified graph")
    scheme.write_buffers[interface] = buffer # Write buffer for message
end

function setWriteBuffer(edge::Edge, buffer::Vector=Array(ProbabilityDistribution,0), scheme::InferenceScheme=currentGraph().active_scheme)
    #(edge in edges(graph)) || error("The specified edge is not part of the current or specified graph")
    scheme.write_buffers[edge] = buffer # Write buffer for marginal
end

function clearBuffers(scheme::InferenceScheme=currentGraph().active_scheme)
    scheme.read_buffers = Dict{TerminalNode, Vector}()
    scheme.write_buffers = Dict{Union(Edge,Interface), Vector}()
    return scheme
end

function setTimeWrap(from::TerminalNode, to::TerminalNode, storage_scheme::InferenceScheme=currentGraph().active_scheme)
    !is(from, to) || error("Cannot create time wrap: from and to must be different nodes")
    push!(storage_scheme.time_wraps, (from, to))
end

clearTimeWraps(scheme::InferenceScheme=graph.currentGraph().active_scheme) = (scheme.time_wraps = Array((TerminalNode, TerminalNode), 0))

function step(scheme::InferenceScheme=currentGraph().active_scheme; n_iterations::Int64=1)
    # Reset marginals
    setVagueMarginals!(scheme)
    # Read buffers
    for (terminal_node, read_buffer) in scheme.read_buffers
        !isempty(read_buffer) || error("Read buffer for node $(terminal_node) is empty")
        terminal_node.value = shift!(read_buffer) # pick the first element off the read_buffer
    end
    # Execute schedule
    for iteration = 1:n_iterations
        execute(scheme)
    end
    # Write buffers
    for (component, write_buffer) in scheme.write_buffers
        if typeof(component) == Interface
            isdefined(component.message, :payload) || error("Cannot write message payload to buffer since there is no message present on $(component)")
            push!(write_buffer, deepcopy(component.message.payload))
        elseif typeof(component) == Edge
            push!(write_buffer, calculateMarginal(component))
        end
    end
    # Time wrap
    for (from, to) in scheme.time_wraps
        isdefined(from.out.partner.message, :payload) || error("There is no message to move to $(to) for the next timestep")
        to.value = deepcopy(from.out.partner.message.payload)
    end
end
step(graph::FactorGraph; n_iterations::Int64=1) = step(graph.active_scheme; n_iterations=n_iterations)