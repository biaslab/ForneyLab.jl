#####################
# Unit tests
#####################

facts("General node properties unit tests") do
    FactorGraph()
    c = 0
    for node_type in [subtypes(Node), subtypes(CompositeNode)]
        if node_type!=CompositeNode && node_type!=MockNode
            context("$(node_type) properties should include interfaces and id") do
                @fact typeof(node_type().interfaces) => Array{Interface, 1} # Check for interface array
                @fact length(node_type().interfaces) >= 1 => true # Check length of interface array
                @fact typeof(node_type().id) => Symbol
                @fact typeof(node_type().i) <: Dict => true
            end

            context("$(node_type) constructor should assign an id") do
                my_node = node_type(;id=symbol("my_id_$(c)"))
                @fact my_node.id => symbol("my_id_$(c)")
            end

            context("$(node_type) constructor should assign interfaces to itself") do
                my_node = node_type()
                for interface_index in 1:length(my_node.interfaces)
                    # Check if the node interfaces couple back to the same node
                    @fact my_node.interfaces[interface_index].node => my_node
                end
            end

            context("$(node_type) should have at least 1 sumProduct!() method") do
                @fact contains(string(methods(ForneyLab.sumProduct!)), string("::", node_type)) => true
            end

            context("$(node_type) constructor should add node to the current graph") do
                my_node = node_type()
                @fact current_graph.n[my_node.id] => my_node
            end

            facts("$(node_type) constructor should check for unique id") do
                MockNode(id=symbol("mock_$(c)"))
                @fact_throws MockNode(id=symbol("mock_$(c)"))
            end
        end
        c += 1
    end
end

#####################
# Integration tests
#####################

facts("Connections between nodes integration tests") do
    context("Nodes can directly be coupled through interfaces by using the interfaces array") do
        initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        n(:node1).interfaces[1].partner = n(:node2).interfaces[1]
        n(:node2).interfaces[1].partner = n(:node1).interfaces[1]
        testInterfaceConnections(n(:node1), n(:node2))
    end

    context("Nodes can directly be coupled through interfaces by using the explicit interface handles") do
        initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        n(:node1).i[:in].partner = n(:node2).i[:out]
        n(:node2).i[:out].partner = n(:node1).i[:in]
        testInterfaceConnections(n(:node1), n(:node2))
    end

    context("Nodes can be coupled by edges by using the interfaces array") do
        initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(n(:node2).interfaces[1], n(:node1).interfaces[1]) # Edge from node 2 to node 1
        testInterfaceConnections(n(:node1), n(:node2))
    end
end
