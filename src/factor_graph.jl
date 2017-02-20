"""
A factor graph consisting of factor nodes and edges.
"""
type FactorGraph
    nodes::Dict{Symbol, Node}
    edges::Dict{Symbol, Edge}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
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
                                            Dict{DataType, Int}())

"""
Automatically generate a unique node id based on the current count of nodes of that type in the graph
"""
function generateNodeId(t::DataType)
    current_graph = currentGraph()
    haskey(current_graph.counters, t) ? current_graph.counters[t] += 1 : current_graph.counters[t] = 1
    count = current_graph.counters[t]
    str = replace(lowercase(split(string(t.name),'.')[end]), "node", "")
    return Symbol("$(str)$(count)")
end

"""
Add a Node to a FactorGraph
"""
function addNode!(graph::FactorGraph, nd::Node)
    !haskey(graph.nodes, nd.id) || error("Graph already contains a Node with id $(nd.id)")
    graph.nodes[nd.id] = nd
    return graph
end