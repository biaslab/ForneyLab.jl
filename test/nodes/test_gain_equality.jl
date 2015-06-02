#####################
# Unit tests
#####################

facts("GainEqualityNode unit tests") do
    context("GainEqualityNode() should initialize a GainEqualityNode with 3 interfaces") do
        FactorGraph()
        GainEqualityNode(id=:node)
        @fact typeof(n(:node)) => GainEqualityNode
        @fact length(n(:node).interfaces) => 3
        @fact n(:node).i[:in1] => n(:node).interfaces[1]
        @fact n(:node).i[:in2] => n(:node).interfaces[2]
        @fact n(:node).i[:out] => n(:node).interfaces[3]
        @fact typeof(n(:node).A) => Array{Float64, 2}
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (xi,W) parametrization") do
        # Forward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:out]))
        @fact n(:gec_node).i[:out].message.payload => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), xi=[1.0, 2.0])
        # Backward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in2]))
        @fact n(:gec_node).i[:in2].message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
        initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), xi=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in1]))
        @fact n(:gec_node).i[:in1].message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), xi=[3.0, 6.0])
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (m,W) parametrization") do
        # Forward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:out]))
        @fact n(:gec_node).i[:out].message.payload => GaussianDistribution(W=reshape([0.5, 0.25, 0.25, 0.5], 2, 2), m=[2.0, 4.0])
        # Backward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in2]))
        @fact n(:gec_node).i[:in2].message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
        initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(W=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in1]))
        @fact n(:gec_node).i[:in1].message.payload => GaussianDistribution(W=reshape([5.0, 2.5, 2.5, 5.0], 2, 2), m=[0.6, 1.2])
    end

    context("GainEqualityNode should propagate a GaussianDistribution with (m,V) parametrization") do
        # Forward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:out]))
        @fact n(:gec_node).i[:out].message.payload => GaussianDistribution(V=reshape([2.0, 1.0, 1.0, 2.0], 2, 2), m=[2.0, 4.0])
        # Backward
        initializeGainEqualityNode(2.0*eye(2), [Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in2]))
        @fact n(:gec_node).i[:in2].message.payload => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
        initializeGainEqualityNode(2.0*eye(2), [nothing, Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0])), Message(GaussianDistribution(V=reshape([1.0, 0.5, 0.5, 1.0], 2, 2), m=[1.0, 2.0]))])
        ForneyLab.execute(ScheduleEntry(n(:gec_node).i[:in1]))
        @fact n(:gec_node).i[:in1].message.payload => GaussianDistribution(V=reshape([0.2, 0.1, 0.1, 0.2], 2, 2), m=[0.6, 1.2])
    end
end