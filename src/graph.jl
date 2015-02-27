export  FactorGraph,
        Subgraph

export  currentGraph,
        setCurrentGraph,
        subgraph,
        setVagueMarginals!,
        clearMessages!,
        nodes,
        edges,
        node,
        factorize!

type Subgraph
    nodes::Set{Node}
    internal_edges::Set{Edge}
    external_edges::Set{Edge}
    internal_schedule::Schedule # Schedule for internal message passing (Dauwels step 2)
    external_schedule::ExternalSchedule # Schedule for updates on nodes connected to external edges (Dauwels step 3)
end

function show(io::IO, subgraph::Subgraph)
    graph = currentGraph()
    println(io, "Subgraph $(findfirst(graph.factorization, subgraph))")
    println(io, " # nodes: $(length(subgraph.nodes))")
    println(io, " # internal edges: $(length(subgraph.internal_edges))")
    println(io, " # external edges: $(length(subgraph.external_edges))")
    println(io, "\nSee also:")
    println(io, " draw(::SubGraph)")
    println(io, " show(nodes(::SubGraph))")
    println(io, " show(edges(::SubGraph))")
end

type FactorGraph
    factorization::Vector{Subgraph} # References to subgraphs
    edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal; also determines the ordering of edges
    approximate_marginals::Dict{(Node, Subgraph), ProbabilityDistribution} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
    time_wraps::Vector{(TerminalNode, TerminalNode)}
end

function show(io::IO, factor_graph::FactorGraph)
    nodes_top = nodes(factor_graph, open_composites=false)
    println(io, "FactorGraph")
    println(io, " # nodes: $(length(nodes_top)) ($(length(nodes(factor_graph))) including child nodes)")
    println(io, " # edges (top level): $(length(edges(nodes_top)))")
    println(io, " # factors: $(length(factor_graph.factorization))")
    println(io, "\nSee also:")
    println(io, " draw(::FactorGraph)")
    println(io, " show(nodes(::FactorGraph))")
    println(io, " show(edges(::FactorGraph))")
end

# Create an empty graph
global current_graph = FactorGraph([Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))],
                                   Dict{Edge, Subgraph}(),
                                   Dict{(Node, Subgraph), ProbabilityDistribution}(),
                                   Dict{TerminalNode, Vector}(),
                                   Dict{Union(Edge,Interface), Vector}(),
                                   Array((TerminalNode, TerminalNode), 0))

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph([Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))],
                                            Dict{Edge, Subgraph}(),
                                            Dict{(Node, Subgraph), ProbabilityDistribution}(),
                                            Dict{TerminalNode, Vector}(),
                                            Dict{Union(Edge,Interface), Vector}(),
                                            Array((TerminalNode, TerminalNode), 0))) # Initialize a new factor graph; automatically sets current_graph

function Subgraph() # Construct and add to current graph
    subgraph = Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
    graph = currentGraph()
    push!(graph.factorization, subgraph)
    return subgraph
end

function ensureMarginal!(node::Node, subgraph::Subgraph, graph::FactorGraph, assign_distribution::DataType)
    # Looks for a marginal in the node-subgraph dictionary.
    # If no marginal is present, it sets and returns a vague distribution.
    # Otherwise, it returns the existing marginal. Used for fast marginal calculations.
    try
        return graph.approximate_marginals[(node, subgraph)]
    catch
        if assign_distribution <: ProbabilityDistribution
            return graph.approximate_marginals[(node, subgraph)] = vague(assign_distribution)
        else
            error("Cannot create a marginal of type $(assign_distribution) since a marginal should be <: ProbabilityDistribution")
        end
    end

end

function conformSubgraph!(subgraph::Subgraph)
    # Updates external edges and nodes field based on internal edges
    subgraph.nodes = Set{Node}()
    subgraph.external_edges = Set{Edge}()

    for internal_edge in subgraph.internal_edges # Collect all nodes in the subgraph
        push!(subgraph.nodes, internal_edge.head.node)
        push!(subgraph.nodes, internal_edge.tail.node)
    end
    subgraph.external_edges = setdiff(edges(subgraph.nodes), subgraph.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges

    return subgraph
end

# Get the subgraph in which internal_edge is internal
subgraph(graph::FactorGraph, internal_edge::Edge) = graph.edge_to_subgraph[internal_edge]
subgraph(internal_edge::Edge) = subgraph(currentGraph(), internal_edge)

function nodesConnectedToExternalEdges(graph::FactorGraph, subgraph::Subgraph)
    # Find nodes connected to external edges
    nodes_connected_to_external = Array(Node, 0)
    for external_edge in subgraph.external_edges
         # Check if head node is in subgraph, connected to external, and not already accounted for
        if (external_edge.head.node in subgraph.nodes) && !(external_edge.head.node in nodes_connected_to_external)
            push!(nodes_connected_to_external, external_edge.head.node)
        end
         # Check if tail node is in subgraph, connected to external, and not already accounted for
        if (external_edge.tail.node in subgraph.nodes) && !(external_edge.tail.node in nodes_connected_to_external)
            push!(nodes_connected_to_external, external_edge.tail.node)
        end
    end
    return nodes_connected_to_external
end
nodesConnectedToExternalEdges(subgraph::Subgraph) = nodesConnectedToExternalEdges(currentGraph(), subgraph)

function setVagueMarginals!(graph::FactorGraph=currentGraph())
    # Sets the vague (almost uninformative) marginals in the graph's approximate marginal dictionary at the appropriate places
    for subgraph in graph.factorization
        external_nodes = nodesConnectedToExternalEdges(graph, subgraph)
        for node in external_nodes
            internal_interfaces = Array(Interface, 0)
            for interface in node.interfaces
                if graph.edge_to_subgraph[interface.edge] == subgraph
                    push!(internal_interfaces, interface)
                end
            end
            # From the distribution types of the marginals on the internal edges we can deduce the unknown marginal types
            if length(internal_interfaces) == 1
                # Univariate
                if internal_interfaces[1].edge.distribution_type != Any
                    marginal_type = internal_interfaces[1].edge.distribution_type
                else
                    error("Unspecified distribution type on edge:\n$(internal_interfaces[1].edge)")
                end
            elseif length(internal_interfaces) == 0
                error("The list of internal edges at node $(node) is empty, check your graph definition.")
            else
                # Multivariate
                internal_incoming_message_types = [interface.edge.distribution_type for interface in internal_interfaces]
                marginal_type = getMarginalType(internal_incoming_message_types...)
            end
            graph.approximate_marginals[(node, subgraph)] = vague(marginal_type)
        end
    end

    return graph
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

function nodes(subgraph::Subgraph; open_composites::Bool=true)
    # Return all nodes in subgraph
    all_nodes = copy(subgraph.nodes)

    if open_composites
        children = Set{Node}()
        for n in all_nodes
            if typeof(n) <: CompositeNode
                union!(children, nodes(n, depth=Inf))
            end
        end
        union!(all_nodes, children)
    end

    return all_nodes
end

function nodes(graph::FactorGraph; open_composites::Bool=true)
    # Return all nodes in graph
    all_nodes = Set{Node}()
    for subgraph in graph.factorization
        union!(all_nodes, subgraph.nodes)
    end

    if open_composites
        children = Set{Node}()
        for n in all_nodes
            if typeof(n) <: CompositeNode
                union!(children, nodes(n, depth=Inf))
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
    edge_set = Set{Edge}()
    for subgraph in graph.factorization
        union!(edge_set, subgraph.internal_edges)
    end
    return edge_set
end
edges(;args...) = edges(currentGraph())

function edges(subgraph::Subgraph; include_external=true)
    if include_external
        return union(subgraph.internal_edges, subgraph.external_edges)
    else
        return subgraph.internal_edges
    end
end

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

function extend(edge_set::Set{Edge})
    # Returns the smallest legal subgraph (connected through deterministic nodes) that includes 'edges'

    edge_cluster = Set{Edge}() # Set to fill with edges in equality cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(edge_cluster, current_edge) # Add to edge cluster
        for node in [current_edge.head.node, current_edge.tail.node] # Check both head and tail node for deterministic type
            if isDeterministic(node)
                for interface in node.interfaces
                    if !is(interface.edge, current_edge) && !(interface.edge in edge_cluster) # Is next level edge not seen yet?
                        push!(edges, interface.edge) # Add to buffer to visit sometime in the future
                    end
                end
            end
        end
    end

    return edge_cluster
end
extend(edge::Edge) = extend(Set{Edge}([edge]))

function factorize!(graph::FactorGraph, edge_set::Set{Edge})
    # The set of internal edges needs to be extended to envelope deterministic nodes
    internal_edges = extend(edge_set)

    # We do not support composite nodes with explicit message passing as node connected to an external edge. All these edges should belong to the same subgraph
    connected_nodes = nodes(internal_edges)
    internal_interfaces = Set{Interface}()
    for edge in internal_edges
        push!(internal_interfaces, edge.head)
        push!(internal_interfaces, edge.tail)
    end

    # Add a subgraph containing the edges specified in internal_edges and conform
    for subgraph in graph.factorization
        setdiff!(subgraph.internal_edges, internal_edges) # Remove edges from existing subgraph
    end
    new_subgraph = Subgraph(copy(connected_nodes), copy(internal_edges), Set{Edge}(), Array(Interface, 0), Array(Node, 0)) # Create subgraph
    push!(graph.factorization, new_subgraph) # Add to current graph
    for internal_edge in internal_edges # Point edges to new subgraph in which they are internal
        graph.edge_to_subgraph[internal_edge] = new_subgraph
    end
    for (subgraph_index, subgraph) in enumerate(graph.factorization)
        # Remove empty subgraphs
        if length(subgraph.internal_edges) == 0
            splice!(graph.factorization, subgraph_index)
            continue
        end
        # Update the external edges and node list
        conformSubgraph!(subgraph)
    end
    return new_subgraph
end
factorize!(internal_edges::Set{Edge}) = factorize!(currentGraph(), internal_edges)
factorize!(internal_edge::Edge) = factorize!(Set{Edge}([internal_edge]))
factorize!(graph::FactorGraph, edge::Edge) = factorize!(graph, Set{Edge}([edge]))
factorize!(graph::FactorGraph, internal_edges::Array{Edge, 1}) = factorize!(graph, Set{Edge}(internal_edges))
factorize!(internal_edges::Array{Edge, 1}) = factorize!(currentGraph(), internal_edges)

function factorize!(graph::FactorGraph)
    # Generate a mean field factorization
    (length(graph.factorization) == 1) || error("Cannot perform mean field factorization on an already factorized graph.")
    edges_to_factor = sort!([e for e in graph.factorization[1].internal_edges]) # Cast to array and sort
    while length(edges_to_factor) > 0 # As long as there are edges to factor
        subgraph = factorize!(graph, edges_to_factor[end])

        # Remove all edges in edge_cluster from edges_to_factor, they have just been added to the same factor
        for e in subgraph.internal_edges
            splice!(edges_to_factor, findfirst(edges_to_factor, e))
        end
    end
    return graph
end
factorize!() = factorize!(currentGraph())
