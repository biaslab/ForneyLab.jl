export  FactorGraph, Wrap

export  currentGraph,
        setCurrentGraph,
        clearMessages!,
        clearWraps,
        wraps,
        nodes,
        edges,
        node,
        edge,
        n,
        e

type FactorGraph
    n::Dict{Symbol, Node} # Nodes
    e::Dict{Symbol, Edge} # Edges
    wraps::Vector{(TerminalNode, TerminalNode)}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
    locked::Bool

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
end

# Create an empty graph
global current_graph = FactorGraph( Dict{Symbol, Node}(),
                                    Dict{Symbol, Edge}(),
                                    Array((TerminalNode, TerminalNode), 0),
                                    Dict{DataType, Int}(),
                                    false,
                                    Dict{TerminalNode, Vector}(),
                                    Dict{Union(Edge,Interface), Vector}())

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph(Dict{Symbol, Node}(),
                                            Dict{Symbol, Edge}(),
                                            Array((TerminalNode, TerminalNode), 0),
                                            Dict{DataType, Int}(),
                                            false,
                                            Dict{TerminalNode, Vector}(),
                                            Dict{Union(Edge,Interface), Vector}())) # Initialize a new factor graph; automatically sets current_graph

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

function generateNodeId(t::DataType)
    # Automatically generates a unique node id based on the current count of nodes of that type in the graph
    haskey(current_graph.counters, t) ? current_graph.counters[t] += 1 : current_graph.counters[t] = 1
    count = current_graph.counters[t]
    str = replace(lowercase(string(t)), "node", "")
    return symbol("$(str)$(count)")
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
node(id::Symbol, c::Int, graph::FactorGraph=currentGraph()) = graph.n[s(id, c)] # Quick concatenated lookup
n = node

edge(id::Symbol, graph::FactorGraph=currentGraph()) = graph.e[id]
edge(id::Symbol, c::Int, graph::FactorGraph=currentGraph()) = graph.e[s(id, c)]
e = edge

wraps(g::FactorGraph=current_graph) = g.wraps

function clearWraps(graph::FactorGraph=current_graph)
    graph.wraps = Array(Wrap, 0)
end

# Wrap type
typealias Wrap (TerminalNode, TerminalNode)
function Wrap(from::TerminalNode, to::TerminalNode, graph::FactorGraph=current_graph)
    !is(from, to) || error("Cannot create wrap: from and to must be different nodes")
    push!(graph.wraps, (from, to))
    return (from, to)
end