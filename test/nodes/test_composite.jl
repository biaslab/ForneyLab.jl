facts("CompositeNode integration tests") do
    (g, t_in, t_out) = initializeCompositeGraph()

    # Wrap graph into CompositeNode
    add_3 = CompositeNode(g, t_in, t_out, id=:add_3)

    context("CompositeNode() should construct a composite node from a FactorGraph") do
        @fact ForneyLab.current_graph.n[:add_3] => add_3
        @fact n(:add_3) => add_3
        @fact add_3.i[:in] => add_3.interfaces[1]
        @fact add_3.i[:out] => add_3.interfaces[2]
        @fact add_3.internal_graph => g
        @fact g.locked => true
        @fact (currentGraph() != g) => true
        @fact add_3.interfaceid_to_terminalnode[1] => t_in
        @fact add_3.interfaceid_to_terminalnode[2] => t_out
        @fact add_3.terminalnode_to_interface[t_in] => add_3.interfaces[1]
        @fact add_3.terminalnode_to_interface[t_out] => add_3.interfaces[2]
        @fact isDeterministic(add_3) => false
        @fact add_3.id => :add_3

        # Connect to higher level graph
        t_in = TerminalNode(5.0) # Use same variable name as a test
        t_out = TerminalNode()
        Edge(t_in, add_3.i[:in])
        Edge(add_3.i[:out], t_out)

        # Verify algorithm execution
        algo = SumProduct.Algorithm(add_3.i[:out])
        @fact run(algo) => Message(DeltaDistribution(8.0))
        t_in.value = 10.0
        @fact run(algo) => Message(DeltaDistribution(13.0))
    end

    context("addRule!() should add a computation rule to a composite node") do
        # Specify a custom rule that performs out = in + 5 instead of out = in + 3
        # A call to sumProduct!(add_3, ...) should yield te result of the custom rule
        internal_in = node(:in, add_3.internal_graph)
        internal_adder = node(:adder, add_3.internal_graph)
        function custom_rule(fields)
            return internal_adder.i[:out].message = Message(DeltaDistribution(mean(internal_in.value)[1]+5.0))
        end
        algo2 = Algorithm(custom_rule)
        addRule!(add_3, add_3.i[:out], sumProduct!, algo2)
        @fact run(algo2) => Message(DeltaDistribution(15.0))
    end

    context("CompositeNode() graph argument should be optional") do
        (g, t_in, t_out) = initializeCompositeGraph()
        add_3_2 = CompositeNode(t_in, t_out, id=:add_3_2)
        @fact typeof(add_3_2) => CompositeNode
    end

    context("CompositeNode should check that all terminal nodes are members of the argument graph") do
        FactorGraph()
        t_odd = TerminalNode()
        Edge(t_odd, TerminalNode())
        (g, t_in, t_out) = initializeCompositeGraph()
        @fact_throws CompositeNode(g, t_odd1, t_out, id=:add_3)
    end
end