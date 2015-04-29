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

externalEdges(sg::Subgraph) = setdiff(edges(nodes(sg, open_composites=false)), sg.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges
nodesConnectedToExternalEdges(sg::Subgraph) = intersect(nodes(externalEdges(sg)), nodes(sg, open_composites=false)) # Nodes connected to external edges are the nodes connected to external edges that are also connected to internal edges

draw(subgraph::Subgraph; args...) = ForneyLab.graphviz(ForneyLab.genDot(nodes(subgraph, open_composites=false), subgraph.internal_edges, external_edges=externalEdges(subgraph)); args...)
drawPdf(subgraph::Subgraph, filename::String) = ForneyLab.dot2pdf(ForneyLab.genDot(nodes(subgraph, open_composites=false), subgraph.internal_edges, external_edges=externalEdges(subgraph)), filename)

function nodes(subgraph::Subgraph; open_composites::Bool=true)
    # Return all nodes in subgraph
    all_nodes = copy(nodes(subgraph.internal_edges))

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

function edges(sg::Subgraph; include_external=true)
    if include_external
        return copy(edges(nodes(sg, open_composites=false)))
    else
        return copy(sg.internal_edges)
    end
end
