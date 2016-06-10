#####################
# Unit tests
#####################

facts("FactorGraph unit tests") do
    context("FactorGraph() should construct a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) --> FactorGraph
        @fact fg.nodes --> Dict{Symbol, Node}()
        @fact fg.edges --> Dict{Symbol, Edge}()
        @fact fg.wraps --> Dict{Symbol, ForneyLab.AbstractWrap}()
        @fact fg.counters --> Dict{DataType, Int}()
        @fact ForneyLab.current_graph --> fg # Global should be set
    end

    context("currentGraph() should return the current graph and create one if necessary") do
        @fact typeof(currentGraph()) --> FactorGraph
        @fact currentGraph() --> ForneyLab.current_graph
        @fact is(currentGraph(), currentGraph()) --> true # Create new graph only if necessary
    end

    context("setCurrentGraph() should set a new current graph object") do
        graph1 = FactorGraph()
        graph2 = FactorGraph()
        @fact currentGraph() --> graph2
        setCurrentGraph(graph1)
        @fact currentGraph() --> graph1
    end

    context("generateNodeId should generate a unique node id") do
        FactorGraph()
        node1 = MockNode()
        node2 = MockNode()
        @fact node1.id --> :mock1
        @fact node2.id --> :mock2
    end
end

#####################
# Integration tests
#####################

import Graphs.num_vertices, Graphs.vertices, Graphs.out_neighbors

facts("FactorGraph related functions") do
    context("edges() should get all edges connected to the node set") do
        initializeLoopyGraph()
        @fact edges(Set{Node}(Node[n(:driver), n(:inhibitor)])) --> Set{Edge}(Edge[n(:driver).i[:in].edge, n(:add).i[:in1].edge, n(:add).i[:out].edge])
    end

    context("nodes() should return a set of all nodes in the graph") do
        # FactorGraph test
        initializeLoopyGraph()
        @fact nodes(currentGraph()) --> Set{Node}([n(:driver), n(:inhibitor), n(:noise), n(:add)])
    end

    context("edge(id::Symbol) should return edge with matching id") do
        graph_1 = FactorGraph()
        graph_2 = FactorGraph()
        my_edge = Edge(MockNode().i[:out], MockNode().i[:out])
        @fact_throws edge(:something)
        @fact_throws edge(:mock1_mock2, graph_1)
        setCurrentGraph(graph_1)
        @fact_throws edge(:mock1_mock2)
        @fact edge(:mock1_mock2, graph_2) --> my_edge
        setCurrentGraph(graph_2)
        @fact edge(:mock1_mock2) --> my_edge
    end

    context("node(id::Symbol) should return node with matching id") do
        graph_1 = FactorGraph()
        graph_2 = FactorGraph()
        mocknode = MockNode(id=:my_mocknode)
        Edge(mocknode.i[:out], MockNode().i[:out])
        @fact_throws node(:some_id)
        @fact_throws node(:my_mocknode, graph_1)
        setCurrentGraph(graph_1)
        @fact_throws node(:my_mocknode)
        @fact node(:my_mocknode, graph_2) --> mocknode
        setCurrentGraph(graph_2)
        @fact node(:my_mocknode) --> mocknode
    end

    context("addNode!() should add a node to a graph") do
        FactorGraph()
        nd = TerminalNode(id=:tnode)
        g = FactorGraph()
        ForneyLab.addNode!(g, nd)
        @fact g.nodes[:tnode] --> nd
    end

    context("summaryDependencyGraph(fg) should build correct dependency graph") do
        initializeLoopyGraph()
        fg = currentGraph()
        dg = summaryDependencyGraph(fg)
        @fact num_vertices(dg) --> 8
        for nd in nodes(fg)
            for iface in nd.interfaces
                @fact iface in vertices(dg) --> true
            end
        end
        @fact length(out_neighbors(n(:noise).i[:out], dg)) --> 0
        @fact length(out_neighbors(n(:add).i[:in1], dg)) --> 2
        @fact n(:noise).i[:out] in out_neighbors(n(:add).i[:in1], dg) --> true
        @fact n(:inhibitor).i[:in] in out_neighbors(n(:add).i[:in1], dg) --> true
    end
end
