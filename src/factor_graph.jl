export  FactorGraph

export  currentGraph,
        setCurrentGraph,
        clearMessages!,
        addNode!,
        hasNode,
        hasEdge,
        nodes,
        edges,
        node,
        edge,
        n,
        e

abstract AbstractWrap

type FactorGraph
    nodes::Dict{Symbol, Node} # Nodes
    edges::Dict{Symbol, Edge} # Edges
    wraps::Dict{Symbol, AbstractWrap}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
    locked::Bool

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union{Edge,Interface}, Vector}
end

function currentGraph()
    # Return currently active FactorGraph.
    # Create one if there is none.
    try
        return current_graph
    catch
        return FactorGraph()
    end
end

setCurrentGraph(graph::FactorGraph) = global current_graph = graph

FactorGraph() = setCurrentGraph(FactorGraph(Dict{Symbol, Node}(),
                                            Dict{Symbol, Edge}(),
                                            Dict{Symbol, AbstractWrap}(),
                                            Dict{DataType, Int}(),
                                            false,
                                            Dict{TerminalNode, Vector}(),
                                            Dict{Union{Edge,Interface}, Vector}())) # Initialize a new factor graph; automatically sets current_graph

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
    current_graph = currentGraph()
    haskey(current_graph.counters, t) ? current_graph.counters[t] += 1 : current_graph.counters[t] = 1
    count = current_graph.counters[t]
    str = replace(lowercase(split(string(t.name),'.')[end]), "node", "")
    return symbol("$(str)$(count)")
end

clearMessages!(graph::FactorGraph = currentGraph()) = map(clearMessages!, nodes(graph))

nodes(graph::FactorGraph = currentGraph()) = Set{Node}(values(graph.nodes))

function nodes(edges::Set{Edge})
    # Return all nodes connected to edges
    connected_nodes = Set{Node}()
    for edge in edges
        push!(connected_nodes, edge.head.node)
        push!(connected_nodes, edge.tail.node)
    end

    return connected_nodes
end

edges(graph::FactorGraph = currentGraph()) = Set{Edge}(values(graph.edges))
edges(node::Node) = Set{Edge}([intf.edge for intf in node.interfaces])
edges(nodeset::Set{Node}) = union(map(edges, nodeset)...)

# Search edge and node by id
node(id::Symbol, graph::FactorGraph = currentGraph()) = graph.nodes[id]
n = node

edge(id::Symbol, graph::FactorGraph = currentGraph()) = graph.edges[id]
e = edge

# Add/remove graph elements
function addNode!(graph::FactorGraph, nd::Node)
    # Add a Node to a FactorGraph
    !graph.locked || error("Cannot add a Node to a locked graph")
    !haskey(graph.nodes, nd.id) || error("Graph already contains a Node with id $(nd.id)")
    graph.nodes[nd.id] = nd

    return graph
end

function Base.delete!(graph::FactorGraph, nd::Node)
    hasNode(graph, nd) || error("Graph does not contain node")
    !graph.locked || error("Cannot delete node from locked graph")

    # Delete wraps
    if typeof(nd) == TerminalNode
        for wr in wraps(nd)
            delete!(graph, wr)
        end
    end

    # Detach read buffers from node
    if haskey(graph.read_buffers, nd)
        detachReadBuffer(nd, graph)
    end

    for iface in nd.interfaces
        # Detach and write buffers from edges/interfaces and delete edges
        if iface.edge != nothing
            delete!(graph, iface.edge)
        end
    end

    # Delete node
    delete!(graph.nodes, nd.id)

    return graph
end

function Base.delete!(graph::FactorGraph, eg::Edge)
    hasEdge(graph, eg) || error("Graph does not contain edge")
    !graph.locked || error("Cannot delete node from locked graph")

    # Decouple buffers
    haskey(graph.write_buffers, eg) && detachWriteBuffer(eg, graph)
    haskey(graph.write_buffers, eg.head) && detachWriteBuffer(eg.head, graph)
    haskey(graph.write_buffers, eg.tail) && detachWriteBuffer(eg.tail, graph)

    # Decouple edge and interfaces
    delete!(graph.edges, eg.id)
    eg.head.partner = nothing
    eg.tail.partner = nothing
    eg.head.edge = nothing
    eg.tail.edge = nothing

    return graph
end


# Check existance of graph elements
hasNode(graph::FactorGraph, nd::Node) = (haskey(graph.nodes, nd.id) && is(graph.nodes[nd.id], nd))
hasEdge(graph::FactorGraph, eg::Edge) = (haskey(graph.edges, eg.id) && is(graph.edges[eg.id], eg))
