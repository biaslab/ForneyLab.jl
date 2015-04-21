module VMP

	# Algorithm
    # # Factorizations for approximate marginals
    # factorization::Vector{Subgraph} # References to the graphs factorization required for anwering this inference question
    # edge_to_subgraph::Dict{Edge, Subgraph} # Fast lookup for edge to subgraph in which edge is internal
    # q_distributions::Dict{Set{Edge}, ProbabilityDistribution} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph

    # internal_edges = graph.edges # Collect internal edges from graph
    # nodes = graph.nodes # Collect nodes from graph 
    # subgraph = Subgraph(internal_edges, Array(Interface, 0)) # Create the first factor
    # edge_to_subgraph = Dict([ie for ie in internal_edges], [subgraph for i=1:length(internal_edges)]) # Mapping for edge to subgraph

end