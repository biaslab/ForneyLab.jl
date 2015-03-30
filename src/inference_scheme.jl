export  InferenceScheme

export  currentScheme,
        setCurrentScheme,
        setVagueMarginals!,
        factorize!,
        subgraph

type InferenceScheme
    # An inference scheme holds all attributes that are required to answer one inference question
    graph::FactorGraph
    factorization::Vector{Subgraph} # References to the graphs factorization required for anwering this inference question
    edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal; also determines the ordering of edges
    approximate_marginals::Dict{(Node, Subgraph), ProbabilityDistribution} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph
    read_buffers::Dict{TerminalNode, Vector}
    write_buffers::Dict{Union(Edge,Interface), Vector}
    time_wraps::Vector{(TerminalNode, TerminalNode)}
end

currentScheme() = current_scheme::InferenceScheme
setCurrentScheme(scheme::InferenceScheme) = global current_scheme = scheme # Set a current_scheme

function InferenceScheme(graph::FactorGraph=currentGraph())
    graph.locked = true # Inference definition has begon; lock graph stucture.
    internal_edges = graph.edges # Collect internal edges from graph
    nodes = graph.nodes # Collect nodes from graph 
    subgraph = Subgraph(nodes, internal_edges, Set{Edge}(), Array(Interface, 0), Array(Node, 0)) # Create the first factor
    edge_to_subgraph = Dict([ie for ie in internal_edges], [subgraph for i=1:length(internal_edges)]) # Mapping for edge to subgraph
    scheme = InferenceScheme(graph,
                             [subgraph],
                             edge_to_subgraph, 
                             Dict{(Node, Subgraph), ProbabilityDistribution}(),
                             Dict{TerminalNode, Vector}(),
                             Dict{Union(Edge,Interface), Vector}(),
                             Array((TerminalNode, TerminalNode), 0))
    setCurrentScheme(scheme)
    return scheme
end


# TODO: review!!
# function show(io::IO, scheme::InferenceScheme)
# 
# end

# Get the subgraph in which internal_edge is internal
subgraph(scheme::InferenceScheme, internal_edge::Edge) = scheme.edge_to_subgraph[internal_edge]
subgraph(internal_edge::Edge) = currentScheme().edge_to_subgraph[internal_edge]

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

function setVagueMarginals!(scheme::InferenceScheme=currentScheme())
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
factorize!(internal_edges::Set{Edge}) = factorize!(currentScheme(), internal_edges)
factorize!(internal_edge::Edge) = factorize!(Set{Edge}([internal_edge]))
factorize!(scheme::InferenceScheme, edge::Edge) = factorize!(scheme, Set{Edge}([edge]))
factorize!(scheme::InferenceScheme, internal_edges::Array{Edge, 1}) = factorize!(scheme, Set{Edge}(internal_edges))
factorize!(internal_edges::Array{Edge, 1}) = factorize!(currentScheme(), internal_edges)

function factorize!(scheme::InferenceScheme=currentScheme())
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