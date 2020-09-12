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

"""
Deterministic edges are associated with a PointMass belief.
Find all edges in the graph that are deterministic.
"""
function deterministicEdges(fg::FactorGraph)
    deterministic_edges = Set{Edge}()
    for node in nodes(fg)
        if isa(node, Clamp)
            union!(deterministic_edges, deterministicEdgeSet(node.i[:out].edge))
        end
    end
    
    return deterministic_edges
end

"""
Extend the root_edge to a deterministically connected set, and determine
which edges in the extended set carry PointMass beliefs.
"""
function deterministicEdgeSet(root_edge::Edge)
    # Collect deterministically connected edges
    edge_set = extend(root_edge, terminate_at_soft_factors=true)
    
    # Collect all interfaces in edge set
    interfaces = Interface[]
    for edge in edge_set
        (edge.a == nothing) || push!(interfaces, edge.a)
        (edge.b == nothing) || push!(interfaces, edge.b)
    end
    
    # Identify a schedule that propagates through the entire edge set
    dg = summaryDependencyGraph(edge_set)
    schedule = children(interfaces, dg)

    # Build dictionary of deterministicness
    is_deterministic = Dict{Interface, Bool}()
    for interface in schedule
        is_deterministic[interface] = isDeterministic(interface, is_deterministic)
    end
    
    # Determine deterministic edges
    deterministic_edges = Set{Edge}()
    for edge in edge_set
        deterministic_a = (edge.a != nothing) && is_deterministic[edge.a]
        deterministic_b = (edge.b != nothing) && is_deterministic[edge.b]

        # An edge is deterministic if any of its interfaces are deterministic
        if deterministic_a || deterministic_b
            push!(deterministic_edges, edge)
        end
    end
    
    return deterministic_edges
end

"""
Determine whether interface sends a PointMass message, based on the connected
node, and the deterministic status of the inbound messages.
"""
function isDeterministic(interface::Interface, is_deterministic::Dict{Interface, Bool})
    # First handle immediate cases
    node = interface.node
    if !isa(node, DeltaFactor)
        # We assume a non-delta factor sends stochastic messages in all directions.
        # This includes the CompositeFactor, for which this is the conservative choice.
        return false
    elseif isa(node, Clamp)
        return true # Clamp nodes only send deterministic messages
    end
    
    # Collect deterministicness of inbounds
    inbounds_deterministic = Bool[]
    for iface in node.interfaces
        if iface != interface
            partner = ultimatePartner(iface)
            push!(inbounds_deterministic, is_deterministic[partner])
        end
    end
    
    # Determine deterministicness of outbound from inbounds
    if isa(node, Equality) && any(inbounds_deterministic)
        return true # Deterministic
    elseif any(.!inbounds_deterministic)
        return false # Delta factor with any stochastic inbounds sends a stochastic message
    elseif all(inbounds_deterministic)
        return true # Delta factor with all inbounds deterministic sends a deterministic message
    else # Fallback
        return false # Stochastic
    end
end