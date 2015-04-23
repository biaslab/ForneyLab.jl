module VMP

using ..ForneyLab

    # Algorithm
    # # Factorizations for approximate marginals
    # factorization::Vector{Subgraph} # References to the graphs factorization required for anwering this inference question
    # edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal
    # q_distributions::Dict{Set{Edge}, ProbabilityDistribution} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph

    # internal_edges = graph.edges # Collect internal edges from graph
    # nodes = graph.nodes # Collect nodes from graph 
    # subgraph = Subgraph(internal_edges, Array(Interface, 0)) # Create the first factor
    # edge_to_subgraph = Dict([ie for ie in internal_edges], [subgraph for i=1:length(internal_edges)]) # Mapping for edge to subgraph

    # In execute:

    # Reset marginals
    #setVagueQDistributions!(algorithm)
    # Execute internal/exernal schedules, ...


    # function nodes(subgraph::Subgraph; open_composites::Bool=true)
    #     # Return all nodes in subgraph
    #     all_nodes = copy(nodes(subgraph.internal_edges))

    #     if open_composites
    #         children = Set{Node}()
    #         for n in all_nodes
    #             if typeof(n) <: CompositeNode
    #                 union!(children, nodes(n, depth=typemax(Int64)))
    #             end
    #         end
    #         union!(all_nodes, children)
    #     end

    #     return all_nodes
    # end

    # function edges(sg::Subgraph; include_external=true)
    #     if include_external
    #         return copy(edges(nodes(sg, open_composites=false)))
    #     else
    #         return copy(sg.internal_edges)
    #     end
    # end

    # draw(subgraph::Subgraph; args...) = graphviz(genDot(subgraph.nodes, subgraph.internal_edges, external_edges=externalEdges(subgraph)); args...)
    # drawPdf(subgraph::Subgraph, filename::String) = dot2pdf(genDot(subgraph.nodes, subgraph.internal_edges, external_edges=externalEdges(subgraph)), filename)


end