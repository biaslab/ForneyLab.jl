export
FactorGraph,
currentGraph,
setCurrentGraph,
nodes,
edges,
Terminal,
ultimatePartner


"""
A factor graph consisting of factor nodes and edges.
"""
mutable struct FactorGraph
    nodes::Dict{Symbol, FactorNode}
    edges::Vector{Edge}
    variables::Dict{Symbol, Variable}
    counters::Dict{String, Int} # Counters for automatic node id assignments
    placeholders::Dict{Clamp, Tuple{Symbol, Int}}
end

"""
Return currently active `FactorGraph`.
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
                                            Dict{Type, Int}(),
                                            Dict{Clamp, Tuple{Symbol, Int}}()))

"""
Automatically generate a unique id based on the current counter value for the element type.
"""
function generateId(t::Type)
    graph = currentGraph()
    tname = lowercase(String(nameof(t))) # Remove module prefix from typename
    counter = haskey(graph.counters, tname) ? graph.counters[tname]+1 : 1
    id = Symbol("$(tname)_$(counter)")

    # Make sure we have a unique id (if we can check it)
    collection =
        if t <: FactorNode
            graph.nodes
        elseif t <: Edge
            graph.edges
        elseif t <: Variable
            graph.variables
        end
    if collection != nothing
        # Increase counter until we have a unique id
        while haskey(collection, id)
            counter += 1
            id = Symbol("$(tname)_$(counter)")
        end
    end

    graph.counters[tname] = counter # save counter

    return id
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
hasNode(graph::FactorGraph, nd::FactorNode) = (haskey(graph.nodes, nd.id) && (graph.nodes[nd.id] === nd))

"""
`hasVariable(graph, var)` checks if `var` is part of `graph`.
"""
hasVariable(graph::FactorGraph, var::Variable) = (haskey(graph.variables, var.id) && (graph.variables[var.id] === var))

nodes(graph::FactorGraph = currentGraph()) = Set{FactorNode}(values(graph.nodes))

function nodes(edgeset::Set{Edge})
    # Return all nodes connected to edgeset
    connected_nodes = Set{FactorNode}()
    for edge in edgeset
        (edge.a == nothing) || push!(connected_nodes, edge.a.node)
        (edge.b == nothing) || push!(connected_nodes, edge.b.node)
    end

    return connected_nodes
end

edges(graph::FactorGraph=currentGraph()) = Set{Edge}(graph.edges)
edges(node::FactorNode) = Set{Edge}([intf.edge for intf in node.interfaces])
edges(nodeset::Set{FactorNode}) = union(Set((edges(node) for node=nodeset))...)

"""
Description:

    Terminal is a special type of node that is only used in the internal
    graph of a CompositeFactor. A Terminal is used to terminate an Edge in the
    internal graph that is linked to an interface of the CompositeFactor.

    A Terminal is linked to an interface of the
    CompositeFactor containing the Terminal.

Interfaces:

    1. out

Construction:

    Terminal(id=:some_id)
"""
mutable struct Terminal <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}
    outer_interface::Interface # Interface of CompositeFactor linked to this Terminal

    function Terminal(out::Variable, outer_interface::Interface; id=generateId(Terminal))
        self = new(id, Array{Interface}(undef, 1), Dict{Symbol,Interface}(), outer_interface)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end


"""
ultimatePartner(interface) finds the 'ultimate partner' of interface.
If interface.partner does not belong to a Terminal, it simply returns
interface.partner. In case of a Terminal node, it finds the first
non-Terminal partner on a higher level factor graph.
"""
function ultimatePartner(interface::Interface)
    if (interface.partner != nothing) && isa(interface.partner.node, Terminal)
        return ultimatePartner(interface.partner.node.outer_interface)
    else
        return interface.partner
    end
end
