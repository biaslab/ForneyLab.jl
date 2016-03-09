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
        eg

abstract AbstractWrap

"""
A factor graph consisting of factor nodes and edges.
"""
type FactorGraph
    nodes::Dict{Symbol, Node} # Nodes
    edges::Dict{Symbol, Edge} # Edges
    wraps::Dict{Symbol, AbstractWrap}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
    locked::Bool
    block_size
    current_section::Int64

    # Connections to the outside world
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union{Edge,Interface}, Vector}

    function FactorGraph(nodes::Dict{Symbol, Node},
                         edges::Dict{Symbol, Edge},
                         wraps::Dict{Symbol, AbstractWrap},
                         counters::Dict{DataType, Int},
                         locked:: Bool,
                         read_buffers::Dict{TerminalNode, Vector},
                         write_buffers::Dict{Union{Edge, Interface}, Vector})

        self = new()
        self.nodes = nodes
        self.edges = edges
        self.wraps = wraps
        self.counters = counters
        self.locked = locked
        self.read_buffers = read_buffers
        self.write_buffers = write_buffers
        self.current_section = 1
        return self
    end
    
end


"""
Return currently active FactorGraph.
Create one if there is none.
"""
function currentGraph()
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
    println(io, " block_size: $(factor_graph.block_size)")
    println(io, " current_section: $(factor_graph.current_section)")
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

edges(graph::FactorGraph=currentGraph()) = Set{Edge}(values(graph.edges))
edges(node::Node) = Set{Edge}([intf.edge for intf in node.interfaces])
edges(nodeset::Set{Node}) = union(map(edges, nodeset)...)

# Search edge and node by id
node(id::Symbol, graph::FactorGraph=currentGraph()) = graph.nodes[id]
nodes(ids::Vector{Symbol}, graph::FactorGraph=currentGraph()) = Node[node(id, graph) for id in ids]
# Shorthands
n(id::Symbol, graph::FactorGraph=currentGraph()) = node(id, graph)
n(ids::Vector{Symbol}, graph::FactorGraph=currentGraph()) = nodes(ids, graph)

edge(id::Symbol, graph::FactorGraph = currentGraph()) = graph.edges[id]
edges(ids::Vector{Symbol}, graph::FactorGraph=currentGraph()) = Edge[edge(id, graph) for id in ids]
# Shorthands
eg(id::Symbol, graph::FactorGraph = currentGraph()) = edge(id, graph)
eg(ids::Vector{Symbol}, graph::FactorGraph=currentGraph()) = edges(ids, graph)

# Add/remove graph elements
function addNode!(graph::FactorGraph, nd::Node)
    # Add a Node to a FactorGraph
    !graph.locked || error("Cannot add a Node to a locked graph")
    !haskey(graph.nodes, nd.id) || error("Graph already contains a Node with id $(nd.id)")
    graph.nodes[nd.id] = nd

    return graph
end

# Check existance of graph elements
hasNode(graph::FactorGraph, nd::Node) = (haskey(graph.nodes, nd.id) && is(graph.nodes[nd.id], nd))
hasEdge(graph::FactorGraph, eg::Edge) = (haskey(graph.edges, eg.id) && is(graph.edges[eg.id], eg))
