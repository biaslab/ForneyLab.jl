export clearMessages!, getNodes, getEdges

# Functions to clear ALL MESSAGES in the graph
clearMessages!(graph::FactorGraph) = map(clearMessages!, getNodes(graph, open_composites=true))
clearMessages!() = clearMessages!(getCurrentGraph())

function addChildNodes!(nodes::Set{Node})
    # Add all child nodes to the nodes set
    composite_nodes_stack = Array(CompositeNode, 0) # Composite nodes to open
    for node in nodes
        if typeof(node) <: CompositeNode
            push!(composite_nodes_stack, node)
        end
    end
    # Open all composite nodes
    while length(composite_nodes_stack) > 0
        composite_node = pop!(composite_nodes_stack)
        for field in names(composite_node)
            if typeof(getfield(composite_node, field)) <: Node
                # Add child
                child_node = getfield(composite_node, field)
                push!(nodes, child_node)
                if typeof(child_node) <: CompositeNode
                    push!(composite_nodes_stack, child_node)
                end
            end
        end
    end
    return nodes
end

function getNodes(subgraph::Subgraph; open_composites::Bool=true)
    all_nodes = copy(subgraph.nodes)

    if open_composites; addChildNodes!(all_nodes); end

    return all_nodes
end

function getNodes(graph::FactorGraph; open_composites::Bool=true)
    all_nodes = Set{Node}()
    for subgraph in graph.factorization
        union!(all_nodes, subgraph.nodes)
    end

    if open_composites; addChildNodes!(all_nodes); end

    return all_nodes
end
getNodes(;args...) = getNodes(getCurrentGraph(); args...)

function getNodes(edges::Set{Edge})
    nodes = Set{Node}()
    for edge in edges
        push!(nodes, edge.head.node)
        push!(nodes, edge.tail.node)
    end
    return nodes
end

function getEdges(graph::FactorGraph)
    # Returns the set of edges in the graph
    edge_set = Set{Edge}()
    for subgraph in graph.factorization
        union!(edge_set, subgraph.internal_edges)
    end
    return edge_set
end
getEdges(;args...) = getEdges(getCurrentGraph())

function getEdges(nodes::Set{Node}; include_external=true)
    # Returns the set of edges connected to nodes, including or excluding external edges
    # An external edge has only head or tail in the interfaces belonging to nodes in the nodes array
    edge_set = Set{Edge}()
    for node in nodes
        for interface in node.interfaces
            if include_external
                if interface.edge!=nothing && ((interface.edge.tail.node in nodes) || (interface.edge.head.node in nodes))
                    push!(edge_set, interface.edge)
                end
            else
                if interface.edge!=nothing && (interface.edge.tail.node in nodes) && (interface.edge.head.node in nodes)
                    push!(edge_set, interface.edge)
                end
            end
        end
    end
    return edge_set
end