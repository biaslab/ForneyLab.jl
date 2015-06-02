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
                for interface_id in 1:length(my_node.interfaces)
                    # Check if the node interfaces couple back to the same node
                    @fact my_node.interfaces[interface_id].node => my_node
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

facts("node(id::Symbol) should return node with matching id") do
    graph_1 = FactorGraph()
    graph_2 = FactorGraph()
    n = MockNode(id=:my_mocknode)
    Edge(n.i[:out], MockNode().i[:out])
    @fact_throws node(:some_id)
    @fact_throws node(:my_mocknode, graph_1)
    setCurrentGraph(graph_1)
    @fact_throws node(:my_mocknode)
    @fact node(:my_mocknode, graph_2) => n
    setCurrentGraph(graph_2)
    @fact node(:my_mocknode) => n
end

facts("n(id::Symbol, c::Int) should return node with matching concatenated id") do
    FactorGraph()
    node = MockNode(id=s(:mock,2))
    @fact n(:mock,2) => node
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

    context("Nodes can directly be coupled through interfaces by using the explicit interface handles") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.i[:in].partner = node2.i[:out]
        node2.i[:out].partner = node1.i[:in]
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can be coupled by edges by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end
end

facts("Functions for collecting nodes") do
    context("nodes() should return a set of all nodes in the graph") do
        # FactorGraph test
        testnodes = initializeLoopyGraph()
        @fact nodes(currentGraph()) => Set{Node}(testnodes)

        # Composite node test
        c_node = CompositeNode(currentGraph())
        @fact nodes(c_node) => Set{Node}(testnodes) # nodes() should return internal nodes of a CompositeNode
    end
end