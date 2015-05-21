#####################
# Unit tests
#####################

facts("GainEqualityNode unit tests") do
    context("GainEqualityNode() should initialize a GainEqualityNode with 3 interfaces") do
        node = GainEqualityNode()
        @fact typeof(node) => GainEqualityNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainEqualityNode() should define an internal Equality and FixedGain node") do
        node = GainEqualityNode([5.0], false)
        @fact typeof(node.equality_node) => EqualityNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainEqualityNode() should point its own interfaces to the internal node interfaces") do
        node = GainEqualityNode([1.0], false)
        @fact node.in1.child => node.equality_node.interfaces[1]
        @fact node.in2.child => node.equality_node.interfaces[3]
        @fact node.out.child => node.fixed_gain_node.out
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (xi,W) parametrization") do
        # Forward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(node.out))
        @fact node.out.message.payload => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), xi=[1.0, 2.0])
        # Backward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in2))
        @fact node.in2.message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
        node = initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in1))
        @fact node.in1.message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (m,W) parametrization") do
        # Forward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(node.out))
        @fact node.out.message.payload => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), m=[2.0, 4.0])
        # Backward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in2))
        @fact node.in2.message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
        node = initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in1))
        @fact node.in1.message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (m,V) parametrization") do
        # Forward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(node.out))
        @fact node.out.message.payload => GaussianDistribution(V=reshape([2.0, 1.0, 1.0, 2.0], 2, 2), m=[2.0, 4.0])
        # Backward
        node = initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in2))
        @fact node.in2.message.payload => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
        node = initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(node.in1))
        @fact node.in1.message.payload => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
    end
end