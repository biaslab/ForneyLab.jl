export Subgraph

type Subgraph
    internal_edges::Set{Edge}
    internal_schedule::Schedule # Schedule for internal message passing
end

function Subgraph()
    Subgraph(Set{Edge}(), Array(Interface, 0))
end

externalEdges(sg::Subgraph) = setdiff(edges(nodes(sg, open_composites=false)), sg.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges
nodesConnectedToExternalEdges(sg::Subgraph) = intersect(nodes(externalEdges(sg)), nodes(sg, open_composites=false)) # Nodes connected to external edges are the nodes connected to external edges that are also connected to internal edges

# Get the subgraph in which internal_edge is internal
subgraph(scheme::DataAwareFactorGraph, internal_edge::Edge) = scheme.edge_to_subgraph[internal_edge]