export FactorGraph, Subgraph
export getCurrentGraph, setCurrentGraph, getSubgraph, setUninformativeMarginals!

# FactorGraph and Subgraph types
type Subgraph
    nodes::Set{Node}
    internal_edges::Set{Edge}
    external_edges::Set{Edge}
    internal_schedule::Schedule # Schedule for internal message passing (Dauwels step 2)
    external_schedule::ExternalSchedule # Schedule for updates on nodes connected to external edges (Dauwels step 3)
end
function show(io::IO, subgraph::Subgraph)
    graph = getCurrentGraph()
    println("Subgraph $(findfirst(graph.factorization, subgraph))")
end

type FactorGraph
    factorization::Array{Subgraph, 1} # References to subgraphs
    edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal; also determines the ordering of edges
    approximate_marginals::Dict{(Node, Subgraph), MessagePayload} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph
end

function show(io::IO, factor_graph::FactorGraph)
    nodes_top = getNodes(factor_graph, open_composites=false)
    println(io, "FactorGraph")
    println(io, " # nodes: $(length(nodes_top)) ($(length(getNodes(factor_graph))) including child nodes)")
    println(io, " # edges (top level): $(length(getEdges(nodes_top)))")
end

# Get and set current graph functions
global current_graph = FactorGraph([Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))], Dict{Edge, Subgraph}(), Dict{(Node, Subgraph), MessagePayload}()) # Initialize a current_graph
getCurrentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph([Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))], Dict{Edge, Subgraph}(), Dict{(Node, Subgraph), MessagePayload}())) # Initialize a new factor graph; automatically sets current_graph
function Subgraph() # Construct and add to current graph
    subgraph = Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
    graph = getCurrentGraph()
    push!(graph.factorization, subgraph)
    return subgraph
end

function getOrCreateMarginal(node::Node, subgraph::Subgraph, graph::FactorGraph, assign_distribution::DataType)
    # Looks for a marginal in the node-subgraph dictionary.
    # When no marginal is present, it sets and returns an uninformative distribution.
    # Otherwise it returns the present approximate marginal. Used for fast marginal calculations.
    if !haskey(graph.approximate_marginals, (node, subgraph))
        if assign_distribution <: ProbabilityDistribution 
            graph.approximate_marginals[(node, subgraph)] = uninformative(assign_distribution)
        else
            error("Unknown assign type argument")
        end
    end
    return graph.approximate_marginals[(node, subgraph)]
end

function conformSubgraph!(subgraph::Subgraph)
    # Updates external edges and nodes field based on internal edges

    subgraph.nodes = Set{Node}()
    subgraph.external_edges = Set{Edge}()

    for internal_edge in subgraph.internal_edges # Collect all nodes in the subgraph
        push!(subgraph.nodes, internal_edge.head.node)
        push!(subgraph.nodes, internal_edge.tail.node)
    end
    subgraph.external_edges = setdiff(getEdges(subgraph.nodes), subgraph.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges

    return subgraph
end

function getSubgraph(graph::FactorGraph, internal_edge::Edge)
    # Returns the subgraph in which internal_edge is internal
    return graph.edge_to_subgraph[internal_edge]
end
getSubgraph(internal_edge::Edge) = getSubgraph(getCurrentGraph(), internal_edge)

function getNodesConnectedToExternalEdges(graph::FactorGraph, subgraph::Subgraph)
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
getNodesConnectedToExternalEdges(subgraph::Subgraph) = getNodesConnectedToExternalEdges(getCurrentGraph(), subgraph)

function setUninformativeMarginals!(graph::FactorGraph=getCurrentGraph())
    # Sets the uninformative marginals in the graph's approximate marginal dictionary at the appropriate places
    for subgraph in graph.factorization
        external_nodes = getNodesConnectedToExternalEdges(graph, subgraph)
        for node in external_nodes
            internal_interfaces = Array(Interface, 0)
            for interface in node.interfaces
                if graph.edge_to_subgraph[interface.edge] == subgraph
                    push!(internal_interfaces, interface)
                end
            end
            # From the invertarization of internal edges and the node type we can deduce the marginal types
            if length(internal_interfaces) == 1
                # Univariate
                marginal_type = getMarginalType(internal_interfaces[1].message_payload_type, internal_interfaces[1].partner.message_payload_type)
            elseif length(internal_interfaces) == 0
                error("The list of internal edges at node $(node) is empty, check your graph definition.")
            else
                # Multivariate
                internal_incoming_message_types = [intf.partner.message_payload_type for intf in internal_interfaces]
                marginal_type = getMarginalType(typeof(node), internal_incoming_message_types...)
            end
            graph.approximate_marginals[(node, subgraph)] = uninformative(marginal_type)
        end
    end

    return graph
end