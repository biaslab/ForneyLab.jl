type RecognitionFactorization
    factors::Array{Subgraph, 1}
    edge_to_subgraph::Dict{Edge, Subgraph}
end
function RecognitionFactorization(graph::FactorGraph=currentGraph())
    # Create an initial subgraph that envelopes the entire graph
    internal_edges = edges(graph)
    sg = Subgraph(internal_edges, Interface[], Node[]) # Build a subgraph that contains all edges in the graph
    edge_to_subgraph = Dict{Edge, Subgraph}()
    for edge in internal_edges # Map all edges to the just created subgraph
        merge!(edge_to_subgraph, Dict{Edge, Subgraph}(edge => sg))
    end
    return RecognitionFactorization([sg], edge_to_subgraph) # Build the new factorization
end

show(io::IO, f::RecognitionFactorization) = println(io, "RecognitionFactorization with $(length(f.factors)) factors")

function initializeVagueRecognitionDistributions(f::RecognitionFactorization, recognition_distribution_types::Dict)
    # Sets the vague (almost uninformative) marginals in the graph's approximate marginal dictionary at the appropriate places

    recognition_distributions = Dict{Tuple{Node, Subgraph}, RecognitionDistribution}()

    for subgraph in f.factors
        (length(subgraph.external_schedule) > 0) || warn("External schedule for subgraph $(subgraph) undefined. Run generateSchedule(...) to generate internal and external schedules.")
        for node in subgraph.external_schedule
            # Collect all internal edges connected to node that belong to subgraph
            internal_edges = Edge[]
            for interface in node.interfaces
                if f.edge_to_subgraph[interface.edge] == subgraph # Edge is internal?
                    push!(internal_edges, interface.edge)
                end
            end

            expanded_recognition_distribution_types = expand(recognition_distribution_types)
            # Look up the recognition distribution type for each cluster
            recognition_type = expanded_recognition_distribution_types[internal_edges[1]]

            # Build recognition_distributions dictionary
            recognition_distributions[(node, subgraph)] = RecognitionDistribution(vague(recognition_type), intersect(edges(node), subgraph.internal_edges))
        end
    end

    return recognition_distributions
end

function resetRecognitionDistributions!(recognition_distributions::Dict{Tuple{Node,Subgraph},RecognitionDistribution}) # Reset
    # Before starting a new iteration, the recognition distributions should be reset to vague
    for recognition_distribution in values(recognition_distributions)
        vague!(recognition_distribution.distribution)
    end
    return recognition_distributions
end

function calculateRecognitionDistribution!(recognition_distributions::Dict{Tuple{Node, Subgraph}, RecognitionDistribution}, node::Node, subgraph::Subgraph, factorization::RecognitionFactorization)
    # Calculate the approximate marginal for node from the perspective of subgraph,
    # and store the result in the scheme.recognition_distributions dictionary.

    recognition_distribution = recognition_distributions[(node, subgraph)]
    if length(recognition_distribution.edges) == 1
        # Update for univariate q
        # When there is only one internal edge, the approximate marginal calculation reduces to the naive marginal update
        internal_edge = first(recognition_distribution.edges) # Extract element
        return calculateRecognitionDistribution!(recognition_distribution, internal_edge.tail.message.payload, internal_edge.head.message.payload)
    end

    # Update for multivariate q
    required_inputs = Array(Any, 0)
    for interface in node.interfaces # Iterate over all edges connected to node
        neighbouring_subgraph = factorization.edge_to_subgraph[interface.edge]
        if neighbouring_subgraph == subgraph # edge is internal
            push!(required_inputs, interface.partner.message)
        else # edge is external
            haskey(recognition_distributions, (node, neighbouring_subgraph)) || error("A required recognition distribution for $(node.id) is not present. Please preset (vague) recognition distributions.")
            push!(required_inputs, recognition_distributions[(node, neighbouring_subgraph)].distribution)
        end
    end
    return calculateRecognitionDistribution!(recognition_distribution, node, required_inputs...)
end

function factorize!(edge_set::Set{Edge}, f::RecognitionFactorization=RecognitionFactorization())
    # Returns a factorization object that specifies a new recognition factorization with edge_set as internal edges

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
factorize!(f::RecognitionFactorization, internal_edge::Edge) = factorize!(Set{Edge}([internal_edge]), f)
factorize!(f::RecognitionFactorization, internal_edges::Vector{Edge}) = factorize!(Set{Edge}(internal_edges), f)

function factorize(recognition_distribution_types::Dict, graph=currentGraph())
    # Encodes factorization of the recognition distribution

    # Single edge example; assigns a Gaussian recognition distribution to edge
    # e1 => GaussianDistribution

    # Column vector of edges example; assigns a Gaussian recognition distribution to each row
    # [e1; e2; e3] => GaussianDistribution

    # Matrix of edges example; assign a MvGaussian distribution to each row
    # [e1 e2 e3; <- cluster
    #  e4 e5 e6;             => MvGaussianDistribution{3}
    #  e7 e8 e9]

    factorization = RecognitionFactorization(graph)
    for key in sort(collect(keys(recognition_distribution_types))) # Iterate over recogition distribution types in a deterministic order
        if typeof(key) <: Edge # Only one edge specified
            factorize!(factorization, key)
        else
            for row in 1:size(key, 1) # Each row in key encodes edges that belong to the same cluster
                factorize!(factorization, vec(key[row, :]))
            end
        end
    end

    return factorization
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
