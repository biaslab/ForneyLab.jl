export factorize!, factorizeMeanField!

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

function factorize!(graph::FactorGraph, edge_set::Set{Edge})
    # The set of internal edges needs to be extended to envelope deterministic nodes
    internal_edges = extend(edge_set)

    # We do not support composite nodes with explicit message passing as node connected to an external edge. All these edges should belong to the same subgraph
    nodes = getNodes(internal_edges)
    internal_interfaces = Set{Interface}()
    for edge in internal_edges
        push!(internal_interfaces, edge.head)
        push!(internal_interfaces, edge.tail)
    end

    # Add a subgraph containing the edges specified in internal_edges and conform
    for subgraph in graph.factorization
        setdiff!(subgraph.internal_edges, internal_edges) # Remove edges from existing subgraph
    end
    new_subgraph = Subgraph(Set{Node}(), copy(internal_edges), Set{Edge}(), Array(Interface, 0), Array(Node, 0)) # Create subgraph
    push!(graph.factorization, new_subgraph) # Add to current graph
    for internal_edge in internal_edges # Point edges to new subgraph in which they are internal
        graph.edge_to_subgraph[internal_edge] = new_subgraph
    end
    for (subgraph_index, subgraph) in enumerate(graph.factorization)
        # Remove empty subgraphs
        if length(subgraph.internal_edges) == 0
            splice!(graph.factorization, subgraph_index)
            continue
        end
        # Update the external edges and node list
        conformSubgraph!(subgraph)
    end
    return new_subgraph
end
factorize!(internal_edges::Set{Edge}) = factorize!(getCurrentGraph(), internal_edges)
factorize!(internal_edge::Edge) = factorize!(Set{Edge}([internal_edge]))
factorize!(internal_edges::Array{Edge, 1}) = factorize!(Set{Edge}(internal_edges))

function factorizeMeanField!(graph::FactorGraph)
    # Generate a mean field factorization
    (length(graph.factorization) == 1) || error("Cannot perform mean field factorization on an already factorized graph.")
    internal_edge_set = copy(graph.factorization[1].internal_edges) # Only one subgraph, so these are all top-level edges in the factor graph
    edges_to_factor = sort([e for e in internal_edge_set]) # Cast to array and sort

    while length(edges_to_factor) > 0 # As long as there are edges to factor
        edge = pop!(edges_to_factor) # Pick an edge to factor
        edge_cluster = extend(edge)
        factorize!(graph, edge_cluster)
        # Remove all edges in edge_cluster from edges_to_factor, they have just been added to the same factor
        edge_not_in_cluster = Bool[!(e in edge_cluster) for e in edges_to_factor]
        edges_to_factor = edges_to_factor[edge_not_in_cluster]
    end
    return graph
end
factorizeMeanField!() = factorizeMeanField!(getCurrentGraph())