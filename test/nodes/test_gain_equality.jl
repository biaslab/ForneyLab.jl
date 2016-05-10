# #####################
# # Unit tests
# #####################

# TODO: non-square A

facts("GainEqualityNode unit tests") do
    context("GainEqualityNode() should initialize a GainEqualityNode with 3 interfaces") do
        FactorGraph()
        GainEqualityNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:in1] --> n(:node).interfaces[1]
        @fact n(:node).i[:in2] --> n(:node).interfaces[2]
        @fact n(:node).i[:out] --> n(:node).interfaces[3]
        @fact typeof(n(:node).gain) --> Array{Float64, 2}
    end

    context("GainEqualityNode provide sumProductRule! for Gaussian (m, V)") do
        # Forward message
        A = [2.0].'
        validateOutboundMessage(GainEqualityNode(A),
                                3,
                                [Message(Gaussian(V=1.0, m=1.0)), Message(Gaussian(V=1.0, m=1.0)), nothing],
                                Gaussian(V=2.0, m=2.0))
        # Backward messages
        validateOutboundMessage(GainEqualityNode(A),
                                2,
                                [Message(Gaussian(V=1.0, m=1.0)), nothing, Message(Gaussian(V=1.0, m=1.0))],
                                Gaussian(V=0.2, m=0.6))
        validateOutboundMessage(GainEqualityNode(A),
                                1,
                                [nothing, Message(Gaussian(V=1.0, m=1.0)), Message(Gaussian(V=1.0, m=1.0))],
                                Gaussian(V=0.2, m=0.6))
    end

    context("GainEqualityNode provide sumProductRule! for MvGaussian with (xi,W) parametrization") do
        # Forward message
        A = 2.0*eye(2)
        validateOutboundMessage(GainEqualityNode(A),
                                3,
                                [Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0])), Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0])), nothing],
                                MvGaussian(W=[0.5 0.25; 0.25 0.5], xi=[1.0, 2.0]))
        # Backward messages
        validateOutboundMessage(GainEqualityNode(A),
                                2,
                                [Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0])), nothing, Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0]))],
                                MvGaussian(W=[5.0 2.5; 2.5 5.0], xi=[3.0, 6.0]))
        validateOutboundMessage(GainEqualityNode(A),
                                1,
                                [nothing, Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0])), Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], xi=[1.0, 2.0]))],
                                MvGaussian(W=[5.0 2.5; 2.5 5.0], xi=[3.0, 6.0]))
    end

    context("GainEqualityNode provide sumProductRule! for MvGaussian with (m,W) parametrization") do
        # Forward message
        A = 2.0*eye(2)
        validateOutboundMessage(GainEqualityNode(A),
                                3,
                                [Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing],
                                MvGaussian(W=[0.5 0.25; 0.25 0.5], m=[2.0, 4.0]))
        # Backward messages
        validateOutboundMessage(GainEqualityNode(A),
                                2,
                                [Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing, Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0]))],
                                MvGaussian(W=[5.0 2.5; 2.5 5.0], m=[0.6, 1.2]))
        validateOutboundMessage(GainEqualityNode(A),
                                1,
                                [nothing, Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(W=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0]))],
                                MvGaussian(W=[5.0 2.5; 2.5 5.0], m=[0.6, 1.2]))
    end

    context("GainEqualityNode provide sumProductRule! for MvGaussian with (m,V) parametrization") do
        # Forward message
        A = 2.0*eye(2)
        validateOutboundMessage(GainEqualityNode(A),
                                3,
                                [Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing],
                                MvGaussian(V=[2.0 1.0; 1.0 2.0], m=[2.0, 4.0]))
        # Backward messages
        validateOutboundMessage(GainEqualityNode(A),
                                2,
                                [Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing, Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0]))],
                                MvGaussian(V=[0.2 0.1; 0.1 0.2], m=[0.6, 1.2]))
        validateOutboundMessage(GainEqualityNode(A),
                                1,
                                [nothing, Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0]))],
                                MvGaussian(V=[0.2 0.1; 0.1 0.2], m=[0.6, 1.2]))
    end

    # context("GainEqualityNode should provide sumProductRule! for non-square A") do
    #     # Forward message
    #     A = [1.0 0.5; -0.5 2.0; 0.5 1.0]
    #     validateOutboundMessage(GainEqualityNode(A),
    #                             3,
    #                             [Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing],
    #                             MvGaussian(...)) # 3D
    #     # Backward messages
    #     validateOutboundMessage(GainEqualityNode(A),
    #                             2,
    #                             [Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), nothing, Message(MvGaussian(V=[1.0 0.5; 0.5 1.0; 2.0 3.0], m=[1.0, 2.0, 3.0]))],
    #                             MvGaussian(...)) # 2D
    #     validateOutboundMessage(GainEqualityNode(A),
    #                             1,
    #                             [nothing, Message(MvGaussian(V=[1.0 0.5; 0.5 1.0], m=[1.0, 2.0])), Message(MvGaussian(V=[1.0 0.5; 0.5 1.0; 2.0 3.0], m=[1.0, 2.0, 3.0]))],
    #                             MvGaussian(...) # 2D
    # end
end
