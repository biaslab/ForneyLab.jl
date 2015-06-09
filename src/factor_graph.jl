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

abstract AbstractWrap

type FactorGraph
    n::Dict{Symbol, Node} # Nodes
    e::Dict{Symbol, Edge} # Edges
    wraps::Dict{Symbol, AbstractWrap}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
    locked::Bool

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
end

# Create an empty graph
global current_graph = FactorGraph( Dict{Symbol, Node}(),
                                    Dict{Symbol, Edge}(),
                                    Dict{Symbol, AbstractWrap}(),
                                    Dict{DataType, Int}(),
                                    false,
                                    Dict{TerminalNode, Vector}(),
                                    Dict{Union(Edge,Interface), Vector}())

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph(Dict{Symbol, Node}(),
                                            Dict{Symbol, Edge}(),
                                            Dict{Symbol, AbstractWrap}(),
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
node(id::Symbol, graph::FactorGraph=current_graph) = graph.n[id]
node(id::Symbol, c::Int, graph::FactorGraph=current_graph) = graph.n[s(id, c)] # Quick concatenated lookup
n = node

edge(id::Symbol, graph::FactorGraph=current_graph) = graph.e[id]
edge(id::Symbol, c::Int, graph::FactorGraph=current_graph) = graph.e[s(id, c)]
e = edge

# Add/remove graph elements
function addNode!(graph::FactorGraph, nd::Node)
    # Add a Node to a FactorGraph
    !graph.locked || error("Cannot add a Node to a locked graph")
    !haskey(graph.n, nd.id) || error("Graph already contains a Node with id $(nd.id)")
    graph.n[nd.id] = nd

    return graph
end

function Base.delete!(graph::FactorGraph, nd::Node)
    hasNode(graph, nd) || error("Graph does not contain node")
    !graph.locked || error("Cannot delete node from locked graph")

    for iface in nd.interfaces
        if iface.edge != nothing
            delete!(graph, iface.edge)
        end
        if haskey(graph.write_buffers, iface)
            detachWriteBuffer(iface, graph)
        end
    end 
    if haskey(graph.read_buffers, nd)
        detachReadBuffer(nd, graph)
    end
    delete!(graph.n, nd.id)
    if typeof(nd) == TerminalNode
        for wr in wraps(nd)
            delete!(graph, wr)
        end
    end

    return graph
end

function Base.delete!(graph::FactorGraph, eg::Edge)
    hasEdge(graph, eg) || error("Graph does not contain edge")
    !graph.locked || error("Cannot delete node from locked graph")

    delete!(graph.e, eg.id)
    if haskey(graph.write_buffers, eg)
        detachWriteBuffer(eg, graph)
    end
    
    return graph
end


# Check existance of graph elements
hasNode(graph::FactorGraph, nd::Node) = (haskey(graph.n, nd.id) && is(graph.n[nd.id], nd))
hasEdge(graph::FactorGraph, eg::Edge) = (haskey(graph.e, eg.id) && is(graph.e[eg.id], eg))
