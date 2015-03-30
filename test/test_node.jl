#####################
# Unit tests
#####################

facts("General node properties unit tests") do
    FactorGraph()
    for node_type in [subtypes(Node), subtypes(CompositeNode)]
        if node_type!=CompositeNode && node_type!=MockNode
            context("$(node_type) properties should include interfaces and name") do
                @fact typeof(node_type().interfaces) => Array{Interface, 1} # Check for interface array
                @fact length(node_type().interfaces) >= 1 => true # Check length of interface array
                @fact typeof(node_type().name) => ASCIIString
            end

            context("$(node_type) constructor should assign a name") do
                my_node = node_type(;name="my_name")
                @fact my_node.name => "my_name"
            end

            context("$(node_type) constructor should assign interfaces to itself") do
                my_node = node_type()
                for interface_id in 1:length(my_node.interfaces)
                    # Check if the node interfaces couple back to the same node
                    @fact my_node.interfaces[interface_id].node => my_node
                end
            end

            context("$(node_type) should have at least 1 sumProduct!() method") do
                @fact contains(string(methods(ForneyLab.sumProduct!)), string("::", node_type)) => true
            end
        end
    end

    for node_type in [subtypes(CompositeNode)]
        context("$(node_type) should have property use_composite_update_rules") do
            @fact node_type().use_composite_update_rules => true || false
        end
    end
end

facts("node(name::String) should return node with matching name") do
    graph_1 = FactorGraph()
    graph_2 = FactorGraph()
    n = MockNode(name="my_mocknode")
    Edge(n.out, MockNode().out)
    @fact_throws node("some_name")
    @fact_throws node("my_mocknode", graph_1)
    setCurrentGraph(graph_1)
    @fact_throws node("my_mocknode")
    @fact node("my_mocknode", graph_2) => n
    setCurrentGraph(graph_2)
    @fact node("my_mocknode") => n
end


#####################
# Integration tests
#####################

facts("Connections between nodes integration tests") do
    context("Nodes can directly be coupled through interfaces by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.interfaces[1].partner = node2.interfaces[1]
        node2.interfaces[1].partner = node1.interfaces[1]
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can directly be coupled through interfaces by using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.in1.partner = node2.out
        node2.out.partner = node1.in1
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can be coupled by edges by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("InferenceScheme should add edge and nodes to current subgraph") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        scheme = InferenceScheme()
        graph = currentGraph()
        @fact edge in scheme.factorization[1].internal_edges => true
        @fact node1 in scheme.factorization[1].nodes => true
        @fact node2 in scheme.factorization[1].nodes => true
    end
end

facts("Functions for collecting nodes") do
    context("nodes() should return an array of all nodes in the graph") do
        # FactorGraph test
        testnodes = initializeLoopyGraph()
        scheme = InferenceScheme()
        found_nodes = nodes(currentGraph())
        @fact length(found_nodes) => length(testnodes)
        for node in testnodes
            @fact node in found_nodes => true
        end

        # Subgraph test
        found_nodes = nodes(scheme.factorization[1])
        @fact length(found_nodes) => length(testnodes)
        for node in testnodes
            @fact node in found_nodes => true
        end

        FactorGraph()
        # Composite node test
        compnode = GainEqualityCompositeNode();
        found_nodes = nodes(compnode, depth=1)
        @fact length(found_nodes) => 2
        @fact compnode.equality_node in found_nodes => true
        @fact compnode.fixed_gain_node in found_nodes => true
    end
end