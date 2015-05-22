import ForneyLab.draw, ForneyLab.drawPdf, ForneyLab.nodes, ForneyLab.edges # Import in order to extend
export draw, drawPdf, nodes, edges

type Subgraph
    internal_edges::Set{Edge}
    internal_schedule::Schedule # Schedule for internal message passing
    external_schedule::Array{Node, 1} # Schedule for marginal updates
end

function Subgraph()
    Subgraph(Set{Edge}(), Array(Interface, 0), Array(Node, 0))
end

function show(io::IO, sg::Subgraph)
    println(io, "Subgraph with $(length(sg.internal_edges)) internal edge(s) and $(length(sg.external_schedule)) node(s) connected to external edges.")
    println(io, "\nSee also:")
    println(io, " draw(::SubGraph)")
    println(io, " show(nodes(::SubGraph))")
    println(io, " show(edges(::SubGraph))")
end

externalEdges(sg::Subgraph) = setdiff(edges(nodes(sg)), sg.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges
nodesConnectedToExternalEdges(sg::Subgraph) = intersect(nodes(externalEdges(sg)), nodes(sg)) # Nodes connected to external edges are the nodes connected to external edges that are also connected to internal edges

draw(subgraph::Subgraph; args...) = ForneyLab.graphviz(ForneyLab.genDot(nodes(subgraph), subgraph.internal_edges, external_edges=externalEdges(subgraph)); args...)
drawPdf(subgraph::Subgraph, filename::String) = ForneyLab.dot2pdf(ForneyLab.genDot(nodes(subgraph), subgraph.internal_edges, external_edges=externalEdges(subgraph)), filename)

nodes(subgraph::Subgraph) = copy(nodes(subgraph.internal_edges))

function edges(sg::Subgraph; include_external=true)
    if include_external
        return copy(edges(nodes(sg)))
    else
        return copy(sg.internal_edges)
    end
end
