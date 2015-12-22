#####################
# Unit tests
#####################

# TODO: non-square A and verify results agains formulas

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
        # Backward messages
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(GaussianDistribution(m=0.0, V=2.0)), Message(GaussianDistribution(m=1.0, V=2.0))],
                                GaussianDistribution(m=0.5, V=1.0))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(GaussianDistribution(m=0.0, V=2.0)), nothing, Message(GaussianDistribution(m=1.0, V=2.0))],
                                GaussianDistribution(m=1.0, V=10.0))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions (m,V)") do
        # Forward message
        A = [2.0 3.0; 3.0 2.0]
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2))), nothing],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))
        # Backward messages
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.7999999999999998, -0.19999999999999996], V=[1.5599999999999996 -1.4399999999999997; -1.4399999999999997 1.5599999999999998]))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2)))],
                                MvGaussianDistribution(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions (m,W)") do
        # Forward message
        A = [2.0 3.0; 3.0 2.0]
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2))), nothing],
                                MvGaussianDistribution(m=[1.0, 2.0], W=[0.35294117647058787 -0.31372549019607887; -0.31372549019607865 0.35294117647058876]))
        # Backward messages
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.8000000000000005, -0.1999999999999995], W=[8.666666666666666 8.0; 8.0 8.666666666666666]))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[1.0, 2.0], W=[0.35294117647058787 -0.31372549019607887; -0.31372549019607865 0.35294117647058876]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions (xi,W)") do
        # Forward message
        A = [2.0 3.0; 3.0 2.0]
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2))), nothing],
                                MvGaussianDistribution(m=[0.5000000000000033, 1.0000000000000024], V=[13.500000000000096 12.000000000000085; 12.000000000000076 13.500000000000062]))
        # Backward messages
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.3999999999999999, -0.09999999999999998], V=[0.7799999999999998 -0.7199999999999999; -0.7199999999999999 0.7799999999999999]))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.5000000000000033, 1.0000000000000024], V=[13.500000000000096 12.000000000000085; 12.000000000000076 13.500000000000062]))
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions (different parametrizations)") do
        # Forward message
        A = [2.0 3.0; 3.0 2.0]
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2))), nothing],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
        # Backward messages
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.3999999999999999, -0.09999999999999998], V=[0.7799999999999998 -0.7199999999999999; -0.7199999999999999 0.7799999999999999]))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2)))],
                                MvGaussianDistribution(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
    end

    # context("GainAdditionNode should provide sumProduct! for non-square A") do
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

#####################
# Integration tests
#####################
