export  FactorGraph,
        InferenceScheme,
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

type InferenceScheme
    # An inference scheme holds all attributes that are required to answer one inference question
    name::ASCIIString
    factorization::Vector{Subgraph} # References to the graphs factorization required for anwering this inference question
    edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal; also determines the ordering of edges
    approximate_marginals::Dict{(Node, Subgraph), ProbabilityDistribution} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
    time_wraps::Vector{(TerminalNode, TerminalNode)}
end

type FactorGraph
    nodes::Set{Node}
    edges::Set{Edge}
    inference_schemes::Array{InferenceScheme, 1}
    active_scheme::InferenceScheme
end

function Subgraph()
    # Construct an empty subgraph.
    # The result is not added to the scheme factorization
    Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
end

function InferenceScheme(;name="unnamed")
    # Construct an empty inference scheme.
    # The result is not added to the graph's inference scheme list
    InferenceScheme(name,
                    [Subgraph()],
                    Dict{Edge, Subgraph}(), 
                    Dict{(Node, Subgraph), ProbabilityDistribution}(),
                    Dict{TerminalNode, Vector}(),
                    Dict{Union(Edge,Interface), Vector}(),
                    Array((TerminalNode, TerminalNode), 0))
end

# Create an empty graph
global current_graph = FactorGraph( Set{Node}(),
                                    Set{Edge}(),
                                    [scheme=InferenceScheme()],
                                    scheme)

currentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph(Set{Node}(),
                                            Set{Edge}(),
                                            [scheme=InferenceScheme()],
                                            scheme)) # Initialize a new factor graph; automatically sets current_graph and active scheme

# activeScheme(graph::FactorGraph=currentGraph()) = graph.active_scheme
# function setActiveScheme(scheme::InferenceScheme, graph::FactorGraph=currentGraph())
#     scheme in graph.inference_schemes || error("Argument inference scheme is not in the argument graph's list of interface schemes") 
#     graph.active_scheme = scheme
# end
# function setActiveScheme(name::ASCIIString, graph::FactorGraph=currentGraph())
#     for scheme in graph.inference_schemes
#         if name == scheme.name
#             graph.active_scheme = scheme
#             return scheme
#         end
#     end
#     error("Inference scheme $(name) not found in argument graph's list of inference schemes")
# end

function InferenceScheme(graph::FactorGraph; name="unnamed")
    # Initializes a new inference scheme with one subgraph, adds it to the graph and sets it as the new active scheme
    scheme = InferenceScheme()

    # Initialize only subgraph from graph definition
    subgraph = scheme.factorization[1]
    subgraph.nodes = nodes(graph)
    subgraph.internal_edges = edges(graph)
    for edge in edges(graph)
        merge!(scheme.edge_to_subgraph, {edge=>subgraph})
    end

    push!(graph.inference_schemes, scheme) # Add scheme to graph inference scheme list
    graph.active_scheme = scheme # Set scheme as the active scheme
    return scheme
end

function Subgraph(scheme::InferenceScheme)
    # Creates a new subgraph and adds it to the inference scheme factorization
    subgraph = Subgraph()
    push!(scheme.factorization, subgraph)
    return subgraph
end

# TODO: review!!
# function show(io::IO, subgraph::Subgraph)
#     graph = currentGraph()
#     println(io, "Subgraph $(findfirst(graph.active_scheme.factorization, subgraph))")
#     println(io, " # nodes: $(length(subgraph.nodes))")
#     println(io, " # internal edges: $(length(subgraph.internal_edges))")
#     println(io, " # external edges: $(length(subgraph.external_edges))")
#     println(io, "\nSee also:")
#     println(io, " draw(::SubGraph)")
#     println(io, " show(nodes(::SubGraph))")
#     println(io, " show(edges(::SubGraph))")
# end



# function show(io::IO, factor_graph::FactorGraph)
#     nodes_top = nodes(factor_graph, open_composites=false)
#     println(io, "FactorGraph")
#     println(io, " # nodes: $(length(nodes_top)) ($(length(nodes(factor_graph))) including child nodes)")
#     println(io, " # edges (top level): $(length(edges(nodes_top)))")
#     println(io, " # inference schemes: $(length(factor_graph.inference_schemes))")
#     println(io, "\nSee also:")
#     println(io, " draw(::FactorGraph)")
#     println(io, " show(nodes(::FactorGraph))")
#     println(io, " show(edges(::FactorGraph))")
# end


function ensureMarginal!(node::Node, subgraph::Subgraph, scheme::InferenceScheme, assign_distribution::DataType)
    # Looks for a marginal in the node-subgraph dictionary.
    # If no marginal is present, it sets and returns a vague distribution.
    # Otherwise, it returns the existing marginal. Used for fast marginal calculations.
    try
        return scheme.approximate_marginals[(node, subgraph)]
    catch
        if assign_distribution <: ProbabilityDistribution
            return scheme.approximate_marginals[(node, subgraph)] = vague(assign_distribution)
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
subgraph(scheme::InferenceScheme, internal_edge::Edge) = scheme.edge_to_subgraph[internal_edge]
subgraph(internal_edge::Edge) = subgraph(currentGraph().active_scheme, internal_edge)

function nodesConnectedToExternalEdges(subgraph::Subgraph)
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

function setVagueMarginals!(scheme::InferenceScheme=currentGraph().active_scheme)
    # Sets the vague (almost uninformative) marginals in the graph's approximate marginal dictionary at the appropriate places
    for subgraph in scheme.factorization
        external_nodes = nodesConnectedToExternalEdges(subgraph)
        for node in external_nodes
            internal_interfaces = Array(Interface, 0)
            for interface in node.interfaces
                if scheme.edge_to_subgraph[interface.edge] == subgraph
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
            scheme.approximate_marginals[(node, subgraph)] = vague(marginal_type)
        end
    end

    return scheme
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
                union!(children, nodes(n, depth=typemax(Int64)))
            end
        end
        union!(all_nodes, children)
    end

    return all_nodes
end

function nodes(graph::FactorGraph; open_composites::Bool=true)
    # Return all nodes in graph
    all_nodes = copy(graph.nodes)

    if open_composites
        children = Set{Node}()
        for n in all_nodes
            if typeof(n) <: CompositeNode
                union!(children, nodes(n, depth=typemax(Int64)))
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
    return copy(graph.edges)
end
edges(;args...) = edges(currentGraph())

function edges(subgraph::Subgraph; include_external=true)
    if include_external
        return union(subgraph.internal_edges, subgraph.external_edges)
    else
        return copy(subgraph.internal_edges)
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

function factorize!(scheme::InferenceScheme, edge_set::Set{Edge})
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
    for subgraph in scheme.factorization
        setdiff!(subgraph.internal_edges, internal_edges) # Remove edges from existing subgraph
    end
    new_subgraph = Subgraph(copy(connected_nodes), copy(internal_edges), Set{Edge}(), Array(Interface, 0), Array(Node, 0)) # Create subgraph
    push!(scheme.factorization, new_subgraph) # Add to current scheme
    for internal_edge in internal_edges # Point edges to new subgraph in which they are internal
        scheme.edge_to_subgraph[internal_edge] = new_subgraph
    end
    for (subgraph_index, subgraph) in enumerate(scheme.factorization)
        # Remove empty subgraphs
        if length(subgraph.internal_edges) == 0
            splice!(scheme.factorization, subgraph_index)
            continue
        end
        # Update the external edges and node list
        conformSubgraph!(subgraph)
    end
    return new_subgraph
end
factorize!(internal_edges::Set{Edge}) = factorize!(currentGraph().active_scheme, internal_edges)
factorize!(internal_edge::Edge) = factorize!(Set{Edge}([internal_edge]))
factorize!(scheme::InferenceScheme, edge::Edge) = factorize!(scheme, Set{Edge}([edge]))
factorize!(scheme::InferenceScheme, internal_edges::Array{Edge, 1}) = factorize!(scheme, Set{Edge}(internal_edges))
factorize!(internal_edges::Array{Edge, 1}) = factorize!(currentGraph().active_scheme, internal_edges)

function factorize!(scheme::InferenceScheme=currentGraph().active_scheme)
    # Generate a mean field factorization
    (length(scheme.factorization) == 1) || error("Cannot perform mean field factorization on an already factorized inference scheme.")
    edges_to_factor = sort!([e for e in scheme.factorization[1].internal_edges]) # Cast to array and sort
    while length(edges_to_factor) > 0 # As long as there are edges to factor
        subgraph = factorize!(scheme, edges_to_factor[end])

        # Remove all edges in edge_cluster from edges_to_factor, they have just been added to the same factor
        for e in subgraph.internal_edges
            splice!(edges_to_factor, findfirst(edges_to_factor, e))
        end
    end
    return scheme
end
factorize!(graph::FactorGraph) = factorize!(graph.active_scheme)
factorize!() = factorize!(currentGraph().active_scheme)
