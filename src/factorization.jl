export factorize!, factorizeMeanField!

function factorize!(graph::FactorGraph, internal_edges::Set{Edge})
    # We do not support composite nodes with explicit message passing as node connected to an external edge. All these edges should belong to the same subgraph
    nodes = getNodes(internal_edges)
    internal_interfaces = Set{Interface}()
    for edge in internal_edges
        push!(internal_interfaces, edge.head)
        push!(internal_interfaces, edge.tail)
    end
    for node in nodes
        if typeof(node) <: CompositeNode && node.use_composite_update_rules == false # Composite node is set to use explicit message passing
            # Check that all interfaces are internal
            for interface in node.interfaces
                if !(interface in internal_interfaces)
                    error("Factorization leads CompositeNode with explicit message passing $(node.name) to become connected to an external edge. Please build the node's internals explicitly and manually define where you wish to place the subgraph borders.")
                end
            end
        end
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
    edges_to_factor = getEdges(graph) # All top-level edges in the factor graph

    while length(edges_to_factor) > 0 # As long as there are edges to factor
        edge = pop!(edges_to_factor) # Pick an edge to factor
        # Check connection to equality node
        if typeof(edge.head.node)==EqualityNode || typeof(edge.tail.node)==EqualityNode
            # Collect all other edges that are connected to this one through equality nodes
            edge_cluster = Set{Edge}() # Set to fill with edges in equality cluster
            connected_edges = Set{Edge}({edge})
            while length(connected_edges) > 0 # As long as there are unchecked edges connected through eq nodes
                current_edge = pop!(connected_edges) # Pick one
                push!(edge_cluster, current_edge) # Add to edge cluster
                for node in [current_edge.head.node, current_edge.tail.node] # Check both head and tail node for EqualityNode
                    if typeof(node) == EqualityNode
                        for interface in node.interfaces
                            if !is(interface.edge, current_edge) && !(interface.edge in edge_cluster) # Is next level edge not seen yet?
                                push!(connected_edges, interface.edge) # Add to buffer to visit sometime in the future
                            end
                        end
                    end
                end
            end
            factorize!(graph, edge_cluster)
            # Remove all edges in edge_cluster from edges_to_factor, they have just been added to the same factor
            setdiff!(edges_to_factor, edge_cluster)
        else
            factorize!(graph, Set{Edge}({edge}))
        end
    end
    return graph
end
factorizeMeanField!() = factorizeMeanField!(getCurrentGraph())