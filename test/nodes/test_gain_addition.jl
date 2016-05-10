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
        @fact typeof(n(:node).gain) --> Array{Float64, 2}
    end

    context("GainAdditionNode should be able to pass Gaussians") do
        # Forward message
        A = 2.0
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(Gaussian(m=0.0, V=1.0)), Message(Gaussian(m=1.0, V=2.0)), nothing],
                                Gaussian(m=1.0, V=6.0))
        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(Gaussian(m=0.0, V=2.0)), Message(Gaussian(m=1.0, V=2.0))],
                                Gaussian(m=0.5, V=1.0))
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(Gaussian(m=0.0, V=2.0)), nothing, Message(Gaussian(m=1.0, V=2.0))],
                                Gaussian(m=1.0, V=10.0))
    end

    context("GainAdditionNode should be able to pass MvGaussians: (m,V) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussian(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing],
                                MvGaussian(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussian(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                MvGaussian(m=[0.8, -0.2], V=[1.56 -1.44; -1.44 1.56]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussian(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                MvGaussian(m=[1.0, 2.0], V=[15.0 12.0; 12.0 15.0]))
    end

    context("GainAdditionNode should be able to pass MvGaussians: (m,W) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussian(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussian(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussian(m=[1.0, 2.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussian(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussian(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[0.8, -0.2], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussian(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussian(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[1.0, 2.0], V=[13.5 12.0; 12.0 13.5]))
    end

    context("GainAdditionNode should be able to pass MvGaussians: (xi,W) parametrization") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussian(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussian(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussian(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussian(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
    end

    context("GainAdditionNode should be able to pass MvGaussians: (different parametrizations)") do
        A = [2.0 3.0; 3.0 2.0]

        # Forward
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                MvGaussian(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))

        # Backward
        validateOutboundMessage(GainAdditionNode(A),
                                1,
                                [nothing, Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78]))

        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussian(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussian(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                MvGaussian(m=[0.5, 1.0], V=[13.5 12.0; 12.0 13.5]))
    end

    # context("GainAdditionNode should provide sumProductRule! for non-square A") do
    #     # Forward message
    #     A = [2.0 3.0; 3.0 2.0; 1.0 2.0]
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             3,
    #                             [Message(MvGaussian(m=[0.0, 0.0], V=eye(2))), Message(MvGaussian(m=[1.0, 2.0, 3.0], V=2.0*eye(3))), nothing],
    #                             MvGaussian()) # 3D
    #     # Backward messages
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             1,
    #                             [nothing, Message(MvGaussian(m=[0.0, 0.0, 1.0], V=eye(3))), Message(MvGaussian(m=[1.0, 2.0], V=2.0*eye(2)))],
    #                             MvGaussian()) # 2D
    #     validateOutboundMessage(GainAdditionNode(A),
    #                             2,
    #                             [Message(MvGaussian(m=[0.0, 0.0], V=eye(2))), nothing, Message(MvGaussian(m=[1.0, 2.0, 3.0], V=2.0*eye(3)))],
    #                             MvGaussian()) # 3D
    # end
end
