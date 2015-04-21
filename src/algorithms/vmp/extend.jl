function extend(edge_set::Set{Edge})
    # Returns the smallest legal subgraph (connected through deterministic nodes) that includes 'edges'

    edge_cluster = Set{Edge}() # Set to fill with edges in equality cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(edge_cluster, current_edge) # Add to edge cluster
        for node in [current_edge.head.node, current_edge.tail.node] # Check both head and tail node for deterministic type
            if isDeterministic(node)
                for interface in node.interfaces
                    if !is(interface.edge, current_edge) && !(interface.edge in edge_cluster) # Is next level edge not seen yet?
                        push!(edges, interface.edge) # Add to buffer to visit sometime in the future
                    end
                end
            end
        end
    end

    return edge_cluster
end
extend(edge::Edge) = extend(Set{Edge}([edge]))