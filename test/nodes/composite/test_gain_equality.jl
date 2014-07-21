#####################
# Unit tests
#####################

facts("GainEqualityCompositeNode unit tests") do
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
end

#####################
# Integration tests
#####################

facts("GainEqualityCompositeNode integration tests") do
    context("Edge can connect a normal node to a GainEqualityCompositeNode") do
        (c_node, node) = initializeTerminalAndGainEqNode()
        Edge(node.out, c_node.in1)
        @fact node.out.partner => c_node.in1 # Set correct partners
        @fact c_node.in1.partner => node.out
        @fact c_node.equality_node.interfaces[1].partner => node.out 
        @fact c_node.in1.child => c_node.equality_node.interfaces[1] # Set child 
    end

    context("A GainEqualityCompositeNode should be able to pass a Gaussian message through its internals") do
        # Forward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), false, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution, GaussianDistribution)
        @fact msg.value.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
        @fact msg.value.xi => [1.0, 2.0]
        # Backward message
        node = initializeGainEqualityCompositeNode(2.0*eye(2), false, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution, GaussianDistribution)
        @fact msg.value.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
        @fact msg.value.xi => [3.0, 6.0]
    end

    context("A GainEqualityCompositeNode should pass a Gaussian message using custom update rules for message passing") do
        # The following tests on the update rules correspond to node 5 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        context("GaussianDistribution with (xi,W) parametrization") do
            # Forward
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution, GaussianDistribution)
            @fact msg.value.W => reshape([0.5, 0.25, 0.25, 0.5], 2, 2)
            @fact msg.value.xi => [1.0, 2.0]
            # Backward message
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution, GaussianDistribution)
            @fact msg.value.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
            @fact msg.value.xi => [3.0, 6.0]
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution, GaussianDistribution)
            @fact msg.value.W => reshape([5.0, 2.5, 2.5, 5.0], 2, 2)
            @fact msg.value.xi => [3.0, 6.0]
        end

        context("GaussianDistribution with (m,W) parametrization") do
            # Forward
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution, GaussianDistribution)
            @fact isApproxEqual(msg.value.W, reshape([0.5, 0.25, 0.25, 0.5], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [2.0, 4.0]) => true
            # Backward message
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution, GaussianDistribution)
            @fact isApproxEqual(msg.value.W, reshape([5.0, 2.5, 2.5, 5.0], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [0.6, 1.2]) => true
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution, GaussianDistribution)
            @fact isApproxEqual(msg.value.W, reshape([5.0, 2.5, 2.5, 5.0], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [0.6, 1.2]) => true
        end

        context("GaussianDistribution with (m,V) parametrization") do
            # Forward
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution, GaussianDistribution)
            ensureMVParametrization!(msg.value)
            @fact isApproxEqual(msg.value.V, reshape([2.0, 1.0, 1.0, 2.0], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [2.0, 4.0]) => true
            # Backward message
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution, GaussianDistribution)
            @fact isApproxEqual(msg.value.V, reshape([0.2, 0.1, 0.1, 0.2], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [0.6, 1.2]) => true
            node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution, GaussianDistribution)
            @fact isApproxEqual(msg.value.V, reshape([0.2, 0.1, 0.1, 0.2], 2, 2)) => true
            @fact isApproxEqual(msg.value.m, [0.6, 1.2]) => true
        end
    end
end