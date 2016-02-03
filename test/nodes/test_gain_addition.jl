#####################
# Unit tests
#####################

# TODO: non-square A

facts("GainAdditionNode unit tests") do
    context("GainAdditionNode() should initialize a GainAdditionNode with 3 interfaces") do
        FactorGraph()
        GainAdditionNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:in1] --> n(:node).interfaces[1]
        @fact n(:node).i[:in2] --> n(:node).interfaces[2]
        @fact n(:node).i[:out] --> n(:node).interfaces[3]
        @fact typeof(n(:node).A) --> Array{Float64, 2}
    end

    context("GainAdditionNode should be able to pass GaussianDistributions") do
        # Forward message
        A = 2.0
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(GaussianDistribution(m=0.0, V=1.0)), Message(GaussianDistribution(m=1.0, V=2.0)), nothing],
                                GaussianDistribution(m=1.0, V=6.0))
        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(GaussianDistribution(m=0.0, V=2.0)), Message(GaussianDistribution(m=1.0, V=2.0))],
                                GaussianDistribution(m=0.5, V=1.0))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(GaussianDistribution(m=0.0, V=2.0)), nothing, Message(GaussianDistribution(m=1.0, V=2.0))],
                                GaussianDistribution(m=1.0, V=10.0))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: (m,V) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.8, -0.2], V=[1.56 -1.44; -1.44 1.56]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: (m,W) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.8, -0.2], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[13.5 12.0; 12.0 13.5]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: (xi,W) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: (different parametrizations)") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
    end

    # context("GainAdditionNode should provide sumProductRule! for non-square A") do
    #     # Forward message
    #     A = [2.0 3.0; 3.0 2.0; 1.0 2.0]
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             3,
    #                             [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), Message(MvGaussianDistribution(m=[1.0, 2.0, 3.0], V=2.0*eye(3))), nothing],
    #                             MvGaussianDistribution()) # 3D
    #     # Backward messages
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             1,
    #                             [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0, 1.0], V=eye(3))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2)))],
    #                             MvGaussianDistribution()) # 2D
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             2,
    #                             [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0, 3.0], V=2.0*eye(3)))],
    #                             MvGaussianDistribution()) # 3D
    # end
end
