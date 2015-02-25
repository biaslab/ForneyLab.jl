#####################
# Unit tests
#####################

facts("General node properties unit tests") do
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

    context("Edge constructor should add edge and nodes to current subgraph") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        graph = currentGraph()
        @fact edge in graph.factorization[1].internal_edges => true
        @fact node1 in graph.factorization[1].nodes => true
        @fact node2 in graph.factorization[1].nodes => true
    end
end