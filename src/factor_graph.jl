export  FactorGraph

export  currentGraph,
        setCurrentGraph,
        clearMessages!,
        nodes,
        edges,
        node

type FactorGraph
    nodes::Set{Node}
    edges::Set{Edge}

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
    wraps::Vector{(TerminalNode, TerminalNode)}
end

# Create an empty graph
global current_graph = FactorGraph( Set{Node}(),
                                    Set{Edge}(),
                                    Dict{TerminalNode, Vector}(),
                                    Dict{Union(Edge,Interface), Vector}(),
                                    Array((TerminalNode, TerminalNode), 0))

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph(Set{Node}(),
                                            Set{Edge}(),
                                            Dict{TerminalNode, Vector}(),
                                            Dict{Union(Edge,Interface), Vector}(),
                                            Array((TerminalNode, TerminalNode), 0))) # Initialize a new factor graph; automatically sets current_graph

function show(io::IO, factor_graph::FactorGraph)
    nodes_top = nodes(factor_graph, open_composites=false)
    println(io, "FactorGraph")
    println(io, " # nodes: $(length(nodes_top)) ($(length(nodes(factor_graph))) including child nodes)")
    println(io, " # edges (top level): $(length(edges(nodes_top)))")
    println(io, "\nSee also:")
    println(io, " draw(::FactorGraph)")
    println(io, " show(nodes(::FactorGraph))")
    println(io, " show(edges(::FactorGraph))")
end

# Functions to clear ALL MESSAGES in the graph
clearMessages!(graph::FactorGraph) = map(clearMessages!, nodes(graph, open_composites=true))
clearMessages!() = clearMessages!(currentGraph())

function nodes(node::CompositeNode; depth::Integer=1)
    # Return set of child nodes up to a certain depth
    # depth = 1 only returns the direct children
    # depth = Inf returns all descendants

    children = Set{Node}()
    composite_nodes_stack = CompositeNode[node] # Composite nodes to open

    generation = 1
    while generation<=depth && length(composite_nodes_stack) > 0
        composite_node = pop!(composite_nodes_stack)
        for field in names(composite_node)
            if typeof(getfield(composite_node, field)) <: Node
                # Add child
                child_node = getfield(composite_node, field)
                push!(children, child_node)
                if typeof(child_node) <: CompositeNode
                    push!(composite_nodes_stack, child_node)
                end
            end
        end
        generation += 1 # keep track of depth in the family tree
    end

    return children
end

function nodes(graph::FactorGraph; open_composites::Bool=true)
    # Return all nodes in graph
    all_nodes = copy(graph.nodes)

    if open_composites
        children = Set{Node}()
        for n in all_nodes
            if typeof(n) <: CompositeNode
                union!(children, nodes(n, depth=typemax(Int64)))
            end
        end
        union!(all_nodes, children)
    end

    return all_nodes
end
nodes(;args...) = nodes(currentGraph(); args...)

function nodes(edges::Set{Edge})
    # Return all nodes connected to edges
    connected_nodes = Set{Node}()
    for edge in edges
        push!(connected_nodes, edge.head.node)
        push!(connected_nodes, edge.tail.node)
    end

    return connected_nodes
end

function edges(graph::FactorGraph)
    # Return the set of all edges in graph
    return copy(graph.edges)
end
edges(;args...) = edges(currentGraph())

edges(node::Node) = Set{Edge}([intf.edge for intf in node.interfaces])

function edges(nodeset::Set{Node}; include_external=true)
    # Return the set of edges connected to nodeset, including or excluding external edges
    # An external edge has only head or tail in the interfaces belonging to nodes in nodeset
    edge_set = Set{Edge}()
    for node in nodeset
        for interface in node.interfaces
            if include_external
                if interface.edge!=nothing && ((interface.edge.tail.node in nodeset) || (interface.edge.head.node in nodeset))
                    push!(edge_set, interface.edge)
                end
            else
                if interface.edge!=nothing && (interface.edge.tail.node in nodeset) && (interface.edge.head.node in nodeset)
                    push!(edge_set, interface.edge)
                end
            end
        end
    end
    return edge_set
end

function node(name::ASCIIString, graph::FactorGraph=currentGraph())
    # Return first node found in graph with same name as argument
    for n in nodes(graph, open_composites=true)
        if n.name == name
            return n
        end
    end

    error("No node with name \"$(name)\" in this FactorGraph")
end