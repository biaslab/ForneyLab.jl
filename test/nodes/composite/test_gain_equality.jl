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

    context("GainEqualityCompositeNode() should allow parameters to be clamped") do
        # Fix in1
        node = GainEqualityCompositeNode(in1=GaussianDistribution())
        @fact typeof(node.in1.partner.node) => ForneyLab.ClampNode
        @fact node.in1.partner.message.value => GaussianDistribution()
        # Fix in2
        node = GainEqualityCompositeNode(in2=GaussianDistribution())
        @fact typeof(node.in2.partner.node) => ForneyLab.ClampNode
        @fact node.in2.partner.message.value => GaussianDistribution()
        # Fix out
        node = GainEqualityCompositeNode(out=GaussianDistribution())
        @fact typeof(node.out.partner.node) => ForneyLab.ClampNode
        @fact node.out.partner.message.value => GaussianDistribution()
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

    context("GainEqualityCompositeNode should propagate a GaussianDistribution with (xi,W) parametrization") do
        # Forward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing])
        ForneyLab.updateNodeMessage!(node.out)
        @fact node.out.message.value => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), xi=[1.0, 2.0])
        # Backward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in2)
        @fact node.in2.message.value => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in1)
        @fact node.in1.message.value => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
    end

    context("GainEqualityCompositeNode should propagate a GaussianDistribution with (m,W) parametrization") do
        # Forward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.updateNodeMessage!(node.out)
        @fact node.out.message.value => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), m=[2.0, 4.0])
        # Backward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in2)
        @fact node.in2.message.value => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in1)
        @fact node.in1.message.value => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
    end

    context("GainEqualityCompositeNode should propagate a GaussianDistribution with (m,V) parametrization") do
        # Forward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.updateNodeMessage!(node.out)
        @fact node.out.message.value => GaussianDistribution(V=reshape([2.0, 1.0, 1.0, 2.0], 2, 2), m=[2.0, 4.0])
        # Backward
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in2)
        @fact node.in2.message.value => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
        node = initializeGainEqualityCompositeNode(2.0*eye(2), true, [nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.updateNodeMessage!(node.in1)
        @fact node.in1.message.value => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
    end
end

#####################
# Integration tests
#####################
