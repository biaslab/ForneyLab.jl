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
        @fact node.in1 => node.equality_node.interfaces[1]
        @fact node.in2 => node.equality_node.interfaces[3]
        @fact node.out => node.fixed_gain_node.out
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
        node = GainEqualityCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
        context("GaussianMessage with (xi,W) parametrization") do
            # Forward message
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])
            inbound_messages[2] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])
            ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            @fact node.out.message.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
            @fact node.out.message.xi => [1.0, 2.0]
            # Backward message
            for outbound_interface = [1,2]
                inbound_messages = Array(GaussianMessage, 3)
                inbound_messages[outbound_interface == 1 ? 2 : 1] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])
                inbound_messages[3] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])
                ForneyLab.updateNodeMessage!(outbound_interface, node, inbound_messages)
                @fact node.interfaces[outbound_interface].message.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
                @fact node.interfaces[outbound_interface].message.xi => [3.0, 6.0]
            end
        end

        context("GaussianMessage with (m,W) parametrization") do
            # Forward message
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])
            inbound_messages[2] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])
            ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            @fact node.out.message.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
            @fact isApproxEqual(node.out.message.m, [2.0, 4.0]) => true
            # Backward message
            for outbound_interface = [1,2]
                inbound_messages = Array(GaussianMessage, 3)
                inbound_messages[outbound_interface == 1 ? 2 : 1] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])
                inbound_messages[3] = GaussianMessage(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])
                ForneyLab.updateNodeMessage!(outbound_interface, node, inbound_messages)
                @fact node.interfaces[outbound_interface].message.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
                @fact node.interfaces[outbound_interface].message.m => [0.6, 1.2]
            end
        end

        context("GaussianMessage with (xi,V) parametrization") do
            # Forward message
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[0.0, 2.0])
            inbound_messages[2] = GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[0.0, 2.0])
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            ensureMVParametrization!(msg)
            @fact isApproxEqual(node.out.message.V, reshape([2.0, 1.0, 1.0, 2.0], 2, 2)) => true
            @fact isApproxEqual(node.out.message.m, [2.0, 4.0]) => true
            # Backward message
            for outbound_interface = [1,2]
                inbound_messages = Array(GaussianMessage, 3)
                inbound_messages[outbound_interface == 1 ? 2 : 1] = GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[0.0, 2.0])
                inbound_messages[3] = GaussianMessage(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])
                msg = ForneyLab.updateNodeMessage!(outbound_interface, node, inbound_messages)
                ensureMVParametrization!(msg)
                @fact isApproxEqual(node.interfaces[outbound_interface].message.V, reshape([0.2, 0.1, 0.1, 0.2], 2, 2)) => true
                @fact isApproxEqual(node.interfaces[outbound_interface].message.m, [0.6, 1.2]) => true
            end
        end
    end
end