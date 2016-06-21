import ForneyLab.draw, ForneyLab.drawPdf, ForneyLab.nodes, ForneyLab.edges # Import in order to extend
export Subgraph
export draw, drawPdf, nodes, edges

type Subgraph
    internal_edges::Set{Edge}
    external_edges::Set{Edge}
    nodes_connected_to_external_edges::Vector{Node} # This is a vector instead of a set, so marginal updates have a fixed order
    internal_schedule::Schedule # Schedule for internal message passing

    Subgraph(internal_edges, external_edges, nodes_connected_to_external_edges) = new(internal_edges, external_edges, nodes_connected_to_external_edges)
end

function show(io::IO, sg::Subgraph)
    println(io, "Subgraph with $(length(sg.internal_edges)) internal edge(s), $(length(sg.external_edges)) external edge(s) and $(length(sg.nodes_connected_to_external_edges)) node(s) connected to external edges.")
    println(io, "\nSee also:")
    println(io, " draw(::SubGraph)")
    println(io, " show(nodes(::SubGraph))")
    println(io, " show(edges(::SubGraph))")
end

draw(subgraph::Subgraph; args...) = ForneyLab.graphviz(ForneyLab.genDot(nodes(subgraph), subgraph.internal_edges, external_edges=subgraph.external_edges); args...)
drawPdf(subgraph::Subgraph, filename::AbstractString) = ForneyLab.dot2pdf(ForneyLab.genDot(nodes(subgraph), subgraph.internal_edges, external_edges=subgraph.external_edges), filename)

nodes(subgraph::Subgraph) = copy(nodes(subgraph.internal_edges)) # Return nodes connected to internal edges

function edges(sg::Subgraph; include_external=true) # Return edges connected to nodes of subgraph
    if include_external
        return copy(edges(nodes(sg)))
    else
        return copy(sg.internal_edges)
    end
end
