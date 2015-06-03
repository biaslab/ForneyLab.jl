#####################
# Unit tests
#####################

facts("FactorGraph unit tests") do
    context("FactorGraph() should initialize a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) => FactorGraph
        @fact fg.n => Dict{Symbol, Node}()
        @fact fg.e => Dict{Symbol, Edge}()
        @fact fg.counters => Dict{DataType, Int}()
        @fact current_graph => fg # Global should be set
    end

    context("currentGraph() should return a current graph object") do
        FactorGraph() # Reset
        @fact currentGraph() => current_graph
        my_graph = currentGraph()  # Make local pointer to global variable
        @fact typeof(my_graph) => FactorGraph
    end

    context("setCurrentGraph() should set a new current graph object") do
        my_first_graph = FactorGraph() # Reset
        my_second_graph = FactorGraph()
        @fact my_first_graph == current_graph => false
        setCurrentGraph(my_first_graph)
        @fact my_first_graph == current_graph => true
    end

    context("generateNodeId should generate a unique node id") do
        FactorGraph()
        node1 = MockNode()
        node2 = MockNode()
        @fact node1.id => :mock1
        @fact node2.id => :mock2
    end
end

#####################
# Integration tests
#####################

facts("Functions for collecting nodes and edges") do
    context("edges() should get all edges connected to the node set") do
        initializeLoopyGraph()
        @fact edges(Set{Node}({n(:driver), n(:inhibitor)})) => Set{Edge}({n(:driver).i[:in].edge, n(:add).i[:in1].edge, n(:add).i[:out].edge})
    end

    context("nodes() should return a set of all nodes in the graph") do
        # FactorGraph test
        initializeLoopyGraph()
        @fact nodes(currentGraph()) => Set{Node}([n(:driver), n(:inhibitor), n(:noise), n(:add)])

        # Composite node test
        c_node = CompositeNode(currentGraph())
        @fact nodes(c_node) => Set{Node}([n(:driver, c_node.internal_graph), n(:inhibitor, c_node.internal_graph), n(:noise, c_node.internal_graph), n(:add, c_node.internal_graph)]) # nodes() should return internal nodes of a CompositeNode
    end

    context("edge(id::Symbol) should return edge with matching id") do
        @fact is(e, edge) => true
        graph_1 = FactorGraph()
        graph_2 = FactorGraph()
        my_edge = Edge(MockNode().i[:out], MockNode().i[:out])
        @fact_throws edge(:something)
        @fact_throws edge(:mock1_mock2, graph_1)
        setCurrentGraph(graph_1)
        @fact_throws edge(:mock1_mock2)
        @fact edge(:mock1_mock2, graph_2) => my_edge
        setCurrentGraph(graph_2)
        @fact edge(:mock1_mock2) => my_edge
    end

    context("node(id::Symbol) should return node with matching id") do
        @fact is(n, node) => true
        graph_1 = FactorGraph()
        graph_2 = FactorGraph()
        mocknode = MockNode(id=:my_mocknode)
        Edge(mocknode.i[:out], MockNode().i[:out])
        @fact_throws node(:some_id)
        @fact_throws node(:my_mocknode, graph_1)
        setCurrentGraph(graph_1)
        @fact_throws node(:my_mocknode)
        @fact node(:my_mocknode, graph_2) => mocknode
        setCurrentGraph(graph_2)
        @fact node(:my_mocknode) => mocknode
    end

end
