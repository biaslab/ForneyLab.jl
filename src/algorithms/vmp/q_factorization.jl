type QFactorization
    factors::Array{Subgraph, 1}
    edge_to_subgraph::Dict{Edge, Subgraph}
end
function QFactorization(graph::FactorGraph=currentGraph())
    # Create an initial subgraph that envelopes the entire graph
    internal_edges = edges(graph)
    sg = Subgraph(internal_edges, Interface[], Node[]) # Build a subgraph that contains all edges in the graph
    edge_to_subgraph = Dict{Edge, Subgraph}()
    for edge in internal_edges # Map all edges to the just created subgraph
        merge!(edge_to_subgraph, Dict{Edge, Subgraph}(edge => sg))
    end
    return QFactorization([sg], edge_to_subgraph) # Build the new factorization
end

show(io::IO, f::QFactorization) = println(io, "QFactorization with $(length(f.factors)) factors")

function vagueQDistributions(f::QFactorization)
    # Sets the vague (almost uninformative) marginals in the graph's approximate marginal dictionary at the appropriate places
    
    q_distributions = Dict{(Node, Subgraph), QDistribution}()

    for subgraph in f.factors
        length(subgraph.external_schedule) > 0 || warn("External schedule for subgraph $(subgraph) undefined. Run generateSchedule(...) to generate internal and external schedules.")
        for node in subgraph.external_schedule
            internal_interfaces = Array(Interface, 0)
            for interface in node.interfaces
                if f.edge_to_subgraph[interface.edge] == subgraph
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
            # Build q_distributions dictionary
            q_distributions[(node, subgraph)] = QDistribution(vague(marginal_type), intersect(edges(node), subgraph.internal_edges))
        end
    end

    return q_distributions
end

function vagueQDistributions!(q_distributions::Dict{Tuple{Node,Subgraph},QDistribution})
    # Before starting a new iteration, the q-distributions should be reset to vague
    for q_distribution in values(q_distributions)
        q_distribution.distribution = vague(typeof(q_distribution.distribution)) # Set to vague
    end
    return q_distributions
end

function calculateQDistribution!(q_distributions::Dict{Tuple{Node, Subgraph}, QDistribution}, node::Node, subgraph::Subgraph, factorization::QFactorization)
    # Calculate the approximate marginal for node from the perspective of subgraph,
    # and store the result in the scheme.q_distributions dictionary.

    q_distribution = q_distributions[(node, subgraph)]
    if length(q_distribution.edges) == 1
        # Update for univariate q
        # When there is only one internal edge, the approximate marginal calculation reduces to the naive marginal update
        internal_edge = first(q_distribution.edges) # Extract element
        return calculateQDistribution!(q_distribution, internal_edge.tail.message.payload, internal_edge.head.message.payload)
    end

    # Update for multivariate q
    required_inputs = Array(Any, 0)
    for interface in node.interfaces # Iterate over all edges connected to node
        neighbouring_subgraph = factorization.edge_to_subgraph[interface.edge]
        if neighbouring_subgraph == subgraph # edge is internal
            push!(required_inputs, interface.partner.message)
        else # edge is external
            haskey(q_distributions, (node, neighbouring_subgraph)) || error("A required q-distribution for $(node.id) is not present. Please preset (vague) q-distributions.")
            push!(required_inputs, q_distributions[(node, neighbouring_subgraph)].distribution)
        end
    end
    return calculateQDistribution!(q_distribution, node, required_inputs...)
end

function factorize!(edge_set::Set{Edge}, f::QFactorization=QFactorization())
    # Returns a factorization object that specifies a new q factorization with edge_set as internal edges

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
    for subgraph in f.factors
        setdiff!(subgraph.internal_edges, internal_edges) # Remove edges from existing subgraphs
    end
    new_subgraph = Subgraph(copy(internal_edges), Interface[], Node[]) # Create subgraph
    push!(f.factors, new_subgraph) # Add to current scheme
    for internal_edge in internal_edges # Point edges to new subgraph in which they are internal
        f.edge_to_subgraph[internal_edge] = new_subgraph
    end
    for (subgraph_index, subgraph) in enumerate(f.factors)
        # Remove empty subgraphs
        if length(subgraph.internal_edges) == 0
            splice!(f.factors, subgraph_index)
            continue
        end
    end
    return f
end
factorize!(internal_edge::Edge, f::QFactorization=QFactorization()) = factorize!(Set{Edge}([internal_edge]), f)
factorize!(internal_edges::Vector{Edge}, f::QFactorization=QFactorization()) = factorize!(Set{Edge}(internal_edges), f)

function factorize(graph::FactorGraph=currentGraph())
    # Generates a mean field factorization based on graph

    f = QFactorization(graph) # Starting point

    edges_to_factor = sort!([e for e in f.factors[1].internal_edges]) # Cast to array and sort
    while !isempty(edges_to_factor) # As long as there are edges to factor
        factorize!(edges_to_factor[end], f)

        # Remove all edges in edge_cluster from edges_to_factor, they have just been added to the same factor
        for e in f.factors[end].internal_edges
            splice!(edges_to_factor, findfirst(edges_to_factor, e))
        end
    end
    return f
end

function extend(edge_set::Set{Edge})
    # Returns the smallest legal subgraph (connected through deterministic nodes) that includes 'edges'

    edge_cluster = Set{Edge}() # Set to fill with edges in equality cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(edge_cluster, current_edge) # Add to edge cluster
        for node in [current_edge.head.node, current_edge.tail.node] # Check both head and tail node for deterministic type
            if ForneyLab.isDeterministic(node)
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
