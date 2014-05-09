facts("GainAdditionCompositeNode") do
    context("GainAdditionCompositeNode() should initialize a GainAdditionCompositeNode with 3 interfaces") do
        node = GainAdditionCompositeNode()
        @fact typeof(node) => GainAdditionCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainAdditionCompositeNode() should define an internal AdditionNode and FixedGainNode") do
        node = GainAdditionCompositeNode([5.0], false)
        @fact typeof(node.addition_node) => AdditionNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainAdditionCompositeNode() should point its own interfaces to the internal node interfaces") do
        node = GainAdditionCompositeNode([1.0], false)
        @fact node.in1 => node.fixed_gain_node.interfaces[1]
        @fact node.in2 => node.addition_node.interfaces[2]
        @fact node.out => node.addition_node.out
    end

    context("GainAdditionCompositeNode should be able to pass GaussianMessages: using shortcut rules or internal graph should yield same result") do
        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->
        #        |_______|

        # TODO: implement tests
    end

    context("A GainAdditionCompositeNode should pass a Gaussian message using custom update rules for message passing") do
        # The following tests on the update rules correspond to node 6 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        node = GainAdditionCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
        context("GaussianMessage with (xi,W) parametrization") do
        end

        context("GaussianMessage with (m,W) parametrization") do
        end

        context("GaussianMessage with (m,V) parametrization") do
        end
    end
end