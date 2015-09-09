#####################
# Unit tests
#####################

facts("General node properties unit tests") do
    FactorGraph()
    c = 0
    for node_type in subtypes(Node)
        if node_type!=CompositeNode && node_type!=MockNode
            context("$(node_type) properties should include interfaces and id") do
                test_node = node_type()
                @fact typeof(test_node.interfaces) => Array{Interface, 1} # Check for interface array
                @fact length(test_node.interfaces) >= 1 => true # Check length of interface array
                @fact typeof(test_node.id) => Symbol
                @fact typeof(test_node.i) <: Dict => true
                @fact_throws deepcopy(test_node)
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
                @fact ForneyLab.current_graph.nodes[my_node.id] => my_node
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

    context("delete! should remove a node and coupled read and write buffers") do
        g = initializeChainOfNodes()
        buff_e = attachWriteBuffer(e(:node1_node2))
        buff_i = attachWriteBuffer(n(:node2).i[:out])
        rd_buff = attachReadBuffer(n(:node1), zeros(3))

        delete!(g, n(:node2))
        @fact haskey(g.nodes, :node2) => false
        @fact haskey(g.edges, :node1_node2) => false
        @fact haskey(g.edges, :node2_node3) => false
        @fact length(g.write_buffers) => 0

        @fact length(g.read_buffers) => 1
        delete!(g, n(:node1))
        @fact haskey(g.nodes, :node1) => false
        @fact length(g.read_buffers) => 0
    end
end

facts("copy(::Node)") do
    g1 = initializePairOfNodes()
    test_edge = Edge(n(:node2).interfaces[1], n(:node1).interfaces[1])
    g2 = FactorGraph() # Add a copy of node2 to a new graph
    node2 = n(:node2, g1)
    node2_copy = copy(node2, id=:node2_copy)
    @fact is(node2, node2_copy) => false
    @fact node2_copy.id => :node2_copy
    # Edges should be removed from copy but not from original
    @fact node2.interfaces[1].edge => test_edge
    @fact node2_copy.interfaces[1].edge => nothing
    @fact node2.interfaces[1].partner => n(:node1, g1).interfaces[1]
    @fact node2_copy.interfaces[1].partner => nothing
    @fact node2_copy in nodes(g1) => false
    @fact node2_copy in nodes(g2) => true
end
