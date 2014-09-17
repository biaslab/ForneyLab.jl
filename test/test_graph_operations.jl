#####################
# Integration tests
#####################

facts("Operations on graph integration tests") do
    context("getNodes() should return an array of all nodes in the graph") do
        nodes = initializeLoopyGraph()
        found_nodes = getNodes(getCurrentGraph())
        @fact length(found_nodes) => length(nodes) # FactorGraph test
        for node in nodes
            @fact node in found_nodes => true
        end

        found_nodes = getNodes(getCurrentGraph().factorization[1]) # Subgraph test
        @fact length(found_nodes) => length(nodes)
        for node in nodes
            @fact node in found_nodes => true
        end
    end

    context("getEdges() should get all edges internal (optionally external as well) to the argument node set") do
        nodes = initializeLoopyGraph()
        @fact getEdges(Set{Node}({nodes[1], nodes[2]}), include_external=false) => Set{Edge}({nodes[1].in1.edge})
        @fact getEdges(Set{Node}({nodes[1], nodes[2]})) => Set{Edge}({nodes[1].in1.edge, nodes[4].in1.edge, nodes[4].out.edge})
    end

    context("addChildNodes!() should add composite node's child nodes to the node array") do
        node = initializeGainEqualityCompositeNode(eye(1), false, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact ForneyLab.addChildNodes!(Set{Node}({node})) => Set{Node}({node, node.equality_node, node.fixed_gain_node})
    end
end