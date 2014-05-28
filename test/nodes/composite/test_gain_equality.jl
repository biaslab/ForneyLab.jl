facts("GainEqualityCompositeNode") do
    context("GainEqualityCompositeNode() should initialize a GainEqualityCompositeNode with 3 interfaces") do
        node = GainEqualityCompositeNode()
        @fact typeof(node) => GainEqualityCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainEqualityCompositeNode() should define an internal Equality and FixedGain node") do
        node = GainEqualityCompositeNode([5.0], false)
        @fact typeof(node.equality_node) => EqualityNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainEqualityCompositeNode() should point its own interfaces to the internal node interfaces") do
        node = GainEqualityCompositeNode([1.0], false)
        @fact node.in1.child => node.equality_node.interfaces[1]
        @fact node.in2.child => node.equality_node.interfaces[3]
        @fact node.out.child => node.fixed_gain_node.out
    end

    context("Edge can connect a normal node to a GainEqualityCompositeNode") do
        # Initialize some nodes
        #             c_node
        #           ------------
        # node      |          |
        # [N]--| |--| |--[=]-| |--|
        #    out in1| in1 |    |
        #           |    ...   |

        c_node = GainEqualityCompositeNode()
        node = ConstantNode()
        edge = Edge(node.out, c_node.in1)
        @fact node.out.partner => c_node.in1 # Set correct partners
        @fact c_node.in1.partner => node.out
        @fact c_node.equality_node.interfaces[1].partner => node.out 
        @fact c_node.in1.child => c_node.equality_node.interfaces[1] # Set child 
    end

    context("A GainEqualityCompositeNode should be able to pass a Gaussian message through its internals") do
        #         _________
        #     in1 |       | in2
        # [N]-----|->[=]<-|-----[N]
        #         |   |   |
        #         |   v   |
        #         |  [A]  |
        #         |___|___|
        #             | out
        #             v

        # Forward message
        node_f = GainEqualityCompositeNode(2.0*eye(2), false) # Flag false to indicate internal message passing
        c_node1_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
        c_node2_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
        Edge(c_node1_f.out, node_f.in1)
        Edge(c_node2_f.out, node_f.in2)
        calculateMessage!(node_f.out)
        @fact node_f.out.message.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
        @fact node_f.out.message.xi => [1.0, 2.0]

        #         _________
        #     in1 |       | in2
        # [N]-----|->[=]<-|-----
        #         |   |   |
        #         |   v   |
        #         |  [A]  |
        #         |___|___|
        #             | out
        #             v
        #            [N]

        # Backward message
        node_b = GainEqualityCompositeNode(2.0*eye(2), false) # Flag false to indicate internal message passing
        c_node1_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
        c_node2_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
        Edge(c_node1_b.out, node_b.in1)
        Edge(node_b.out, c_node2_b.out)
        calculateMessage!(node_b.in2)
        @fact node_b.in2.message.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
        @fact node_b.in2.message.xi => [3.0, 6.0]
    end

    context("A GainEqualityCompositeNode should pass a Gaussian message using custom update rules for message passing") do
        # The following tests on the update rules correspond to node 5 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        context("GaussianMessage with (xi,W) parametrization") do
            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v

            # Forward message
            node_f = GainEqualityCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
            c_node1_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            c_node2_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            Edge(c_node1_f.out, node_f.in1)
            Edge(c_node2_f.out, node_f.in2)
            calculateMessage!(node_f.out)
            @fact node_f.out.message.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
            @fact node_f.out.message.xi => [1.0, 2.0]

            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in1)
            Edge(node_b.out, c_node2_b.out)
            calculateMessage!(node_b.in2)
            @fact node_b.in2.message.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
            @fact node_b.in2.message.xi => [3.0, 6.0]

            #         _________
            #     in1 |       | in2
            #    -----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in2)
            Edge(node_b.out, c_node2_b.out)
            calculateMessage!(node_b.in1)
            @fact node_b.in1.message.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
            @fact node_b.in1.message.xi => [3.0, 6.0]
        end

        context("GaussianMessage with (m,W) parametrization") do
            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v

            # Forward message
            node_f = GainEqualityCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
            c_node1_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_f = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_f.out, node_f.in1)
            Edge(c_node2_f.out, node_f.in2)
            calculateMessage!(node_f.out)
            @fact isApproxEqual(node_f.out.message.W, reshape([0.5, 0.25, 0.25, 0.5], 2, 2)) => true
            @fact isApproxEqual(node_f.out.message.m, [2.0, 4.0]) => true

            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in1)
            Edge(node_b.out, c_node2_b.out)
            calculateMessage!(node_b.in2)
            @fact isApproxEqual(node_b.in2.message.W, reshape([5.0, 2.5, 2.5, 5.0], 2, 2)) => true
            @fact isApproxEqual(node_b.in2.message.m, [0.6, 1.2]) => true

            #         _________
            #     in1 |       | in2
            #    -----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in2)
            Edge(node_b.out, c_node2_b.out)
            calculateMessage!(node_b.in1)
            @fact isApproxEqual(node_b.in1.message.W, reshape([5.0, 2.5, 2.5, 5.0], 2, 2)) => true
            @fact isApproxEqual(node_b.in1.message.m, [0.6, 1.2]) => true
        end

        context("GaussianMessage with (xi,V) parametrization") do
            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v

            # Forward message
            node_f = GainEqualityCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
            c_node1_f = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_f = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_f.out, node_f.in1)
            Edge(c_node2_f.out, node_f.in2)
            ensureMVParametrization!(calculateMessage!(node_f.out))
            @fact isApproxEqual(node_f.out.message.V, reshape([2.0, 1.0, 1.0, 2.0], 2, 2)) => true
            @fact isApproxEqual(node_f.out.message.m, [2.0, 4.0]) => true

            #         _________
            #     in1 |       | in2
            # [N]-----|->[=]<-|-----
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in1)
            Edge(node_b.out, c_node2_b.out)
            ensureMVParametrization!(calculateMessage!(node_b.in2))
            @fact isApproxEqual(node_b.in2.message.V, reshape([0.2, 0.1, 0.1, 0.2], 2, 2)) => true
            @fact isApproxEqual(node_b.in2.message.m, [0.6, 1.2]) => true

            #         _________
            #     in1 |       | in2
            #    -----|->[=]<-|-----[N]
            #         |   |   |
            #         |   v   |
            #         |  [A]  |
            #         |___|___|
            #             | out
            #             v
            #            [N]

            # Backward message
            node_b = GainEqualityCompositeNode(2.0*eye(2))
            c_node1_b = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            c_node2_b = ConstantNode(GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))
            Edge(c_node1_b.out, node_b.in2)
            Edge(node_b.out, c_node2_b.out)
            ensureMVParametrization!(calculateMessage!(node_b.in1))
            @fact isApproxEqual(node_b.in1.message.V, reshape([0.2, 0.1, 0.1, 0.2], 2, 2)) => true
            @fact isApproxEqual(node_b.in1.message.m, [0.6, 1.2]) => true
        end
    end
end