#####################
# Unit tests
#####################

facts("Subgraph() should initialize a subgraph on which nodes() and edges() are defined") do
    #               [T2]   [T3]
    #     here       |      | 
    #      v         v      v 
    # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

    g = initializeFactoringGraph()

    internal_edges = Set([eg(:t1_a1), eg(:a1_g1)])
    external_edges = Set([eg(:t2_g1), eg(:g1_g2)])
    nodes_connected_to_external_edges = [n(:g1)]
    sg = ForneyLab.Subgraph(internal_edges, external_edges, nodes_connected_to_external_edges)

    # Test subgraph construction
    @fact typeof(sg) --> ForneyLab.Subgraph
    @fact sg.internal_edges --> internal_edges
    @fact sg.external_edges --> external_edges
    @fact sg.nodes_connected_to_external_edges --> nodes_connected_to_external_edges

    # Test nodes and edges
    @fact nodes(sg) --> Set([n(:t1), n(:a1), n(:g1)])
    @fact edges(sg) --> Set([eg(:t1_a1), eg(:a1_g1), eg(:t2_g1), eg(:g1_g2)])
    @fact edges(sg, include_external=false) --> Set([eg(:t1_a1), eg(:a1_g1)])
end
