export
FactorGraph,
currentGraph,
Terminal


"""
A factor graph consisting of factor nodes and edges.
"""
type FactorGraph
    nodes::Dict{Symbol, FactorNode}
    edges::Vector{Edge}
    variables::Dict{Symbol, Variable}
    counters::Dict{DataType, Int} # Counters for automatic node id assignments
    placeholders::Dict{Constant, Tuple{Symbol, Int}}
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

FactorGraph() = setCurrentGraph(FactorGraph(Dict{Symbol, FactorNode}(),
                                            Edge[],
                                            Dict{Symbol, Variable}(),
                                            Dict{DataType, Int}(),
                                            Dict{Constant, Tuple{Symbol, Int}}()))

"""
Automatically generate a unique id based on the current counter value for the element type.
"""
function generateId(t::DataType)
    current_graph = currentGraph()
    haskey(current_graph.counters, t) ? current_graph.counters[t] += 1 : current_graph.counters[t] = 1
    count = current_graph.counters[t]
    str = lowercase(split(string(t.name),'.')[end]) # Remove "ForneyLab." from typename
    return Symbol("$(str)_$(count)")
end

"""
Add a FactorNode to a FactorGraph
"""
function addNode!(graph::FactorGraph, nd::FactorNode)
    !haskey(graph.nodes, nd.id) || error("Graph already contains a FactorNode with id $(nd.id)")
    graph.nodes[nd.id] = nd
    return graph
end

"""
Add a Variable to a FactorGraph
"""
function addVariable!(graph::FactorGraph, var::Variable)
    !haskey(graph.variables, var.id) || error("Graph already contains a Variable with id $(var.id)")
    graph.variables[var.id] = var
    return graph
end

"""
`hasNode(graph, node)` checks if `node` is part of `graph`.
"""
hasNode(graph::FactorGraph, nd::FactorNode) = (haskey(graph.nodes, nd.id) && is(graph.nodes[nd.id], nd))

"""
`hasVariable(graph, var)` checks if `var` is part of `graph`.
"""
hasVariable(graph::FactorGraph, var::Variable) = (haskey(graph.variables, var.id) && is(graph.variables[var.id], var))

"""
Description:

    Terminal is a special node to terminate an Edge.

Interfaces:

    1. out

Construction:

    Terminal(id=:some_id)
"""
type Terminal <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Terminal(out::Variable; id=generateId(Terminal))
        self = new(id, Array(Interface, 1), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end