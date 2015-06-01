export  FactorGraph

export  currentGraph,
        setCurrentGraph,
        clearMessages!,
        nodes,
        edges,
        node,
        edge,
        n,
        e

type FactorGraph
    n::Dict{Symbol, Node} # Nodes
    e::Dict{Symbol, Edge} # Edges
    locked::Bool

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
    wraps::Vector{(TerminalNode, TerminalNode)}
end

# Create an empty graph
global current_graph = FactorGraph( Dict{Symbol, Node}(),
                                    Dict{Symbol, Edge}(),
                                    false,
                                    Dict{TerminalNode, Vector}(),
                                    Dict{Union(Edge,Interface), Vector}(),
                                    Array((TerminalNode, TerminalNode), 0))

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph(Dict{Symbol, Node}(),
                                            Dict{Symbol, Edge}(),
                                            false,
                                            Dict{TerminalNode, Vector}(),
                                            Dict{Union(Edge,Interface), Vector}(),
                                            Array((TerminalNode, TerminalNode), 0))) # Initialize a new factor graph; automatically sets current_graph

function show(io::IO, factor_graph::FactorGraph)
    println(io, "FactorGraph")
    println(io, " # nodes: $(length(nodes(factor_graph)))")
    println(io, " # edges: $(length(edges(factor_graph)))")
    println(io, " # wraps: $(length(wraps(factor_graph)))")
    println(io, "\nSee also:")
    println(io, " draw(::FactorGraph)")
    println(io, " show(nodes(::FactorGraph))")
    println(io, " show(edges(::FactorGraph))")
    println(io, " show(wraps(::FactorGraph))")
end

clearMessages!(graph::FactorGraph = current_graph) = map(clearMessages!, nodes(graph))

nodes(graph::FactorGraph = current_graph) = Set{Node}(values(graph.n))

function nodes(edges::Set{Edge})
    # Return all nodes connected to edges
    connected_nodes = Set{Node}()
    for edge in edges
        push!(connected_nodes, edge.head.node)
        push!(connected_nodes, edge.tail.node)
    end

    return connected_nodes
end

edges(graph::FactorGraph = current_graph) = Set{Edge}(values(graph.e))
edges(node::Node) = Set{Edge}([intf.edge for intf in node.interfaces])
edges(nodeset::Set{Node}) = union(map(edges, nodeset)...)

# Search edge and node by id
node(id::Symbol, graph::FactorGraph=currentGraph()) = graph.n[id]
n(id::Symbol, graph::FactorGraph=currentGraph()) = node(id, graph)
edge(id::Symbol, graph::FactorGraph=currentGraph()) = graph.e[id]
e(id::Symbol, graph::FactorGraph=currentGraph()) = edge(id, graph)
