#####################
# Unit tests
#####################

facts("AdditionNode unit tests") do
    context("AdditionNode() should initialize an AdditionNode with 3 interfaces") do
        FactorGraph()
        AdditionNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:in1] --> n(:node).interfaces[1]
        @fact n(:node).i[:in2] --> n(:node).interfaces[2]
        @fact n(:node).i[:out] --> n(:node).interfaces[3]
    end

    context("AdditionNode should provide sumProductRule! for Delta{Float64}") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(Delta(2.0)), Message(Delta(3.0)), nothing],
                                Delta(5.0))
        # Backward message
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(Delta(2.0)), Message(Delta(3.0))],
                                Delta(1.0))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(Delta(2.0)), nothing, Message(Delta(3.0))],
                                Delta(1.0))
    end

    context("AdditionNode should provide sumProductRule! for MvDelta{Float64}") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(MvDelta([1.0, 2.0])), Message(MvDelta([3.0, 4.0])), nothing],
                                MvDelta([4.0, 6.0]))
        # Backward message
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(MvDelta([1.0, 2.0])), Message(MvDelta([3.0, 4.0]))],
                                MvDelta([2.0, 2.0]))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(MvDelta([1.0, 2.0])), nothing, Message(MvDelta([3.0, 4.0]))],
                                MvDelta([2.0, 2.0]))
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should provide sumProductRule! for Gaussian") do
        context("Gaussian with (m,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(Gaussian(m=1.0, V=2.0)), Message(Gaussian(m=3.0, V=4.0)), nothing],
                                    Gaussian(m=4.0, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=1.0, V=2.0)), Message(Gaussian(m=3.0, V=4.0))],
                                    Gaussian(m=2.0, V=6.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(Gaussian(m=1.0, V=2.0)), nothing, Message(Gaussian(m=3.0, V=4.0))],
                                    Gaussian(m=2.0, V=6.0))
        end

        context("Gaussian with (m,W) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(Gaussian(m=1.0, W=2.0)), Message(Gaussian(m=3.0, W=4.0)), nothing],
                                    Gaussian(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=1.0, W=2.0)), Message(Gaussian(m=3.0, W=4.0))],
                                    Gaussian(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(Gaussian(m=1.0, W=2.0)), nothing, Message(Gaussian(m=3.0, W=4.0))],
                                    Gaussian(m=2.0, W=4.0/3.0))
        end

        context("Gaussian with (xi,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(Gaussian(xi=1.0, V=2.0)), Message(Gaussian(xi=3.0, V=4.0)), nothing],
                                    Gaussian(xi=14/6, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(Gaussian(xi=1.0, V=2.0)), Message(Gaussian(xi=3.0, V=4.0))],
                                    Gaussian(xi=10/6, V=6.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(Gaussian(xi=1.0, V=2.0)), nothing, Message(Gaussian(xi=3.0, V=4.0))],
                                    Gaussian(xi=10/6, V=6.0))
        end

        context("Gaussian with different parametrizations") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(Gaussian(m=1.0, V=0.5)), Message(Gaussian(m=3.0, W=4.0)), nothing],
                                    Gaussian(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=1.0, V=0.5)), Message(Gaussian(m=3.0, W=4.0))],
                                    Gaussian(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(Gaussian(m=1.0, V=0.5)), nothing, Message(Gaussian(m=3.0, W=4.0))],
                                    Gaussian(m=2.0, W=4.0/3.0))
        end
        context("Support for improper Gaussians") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(Gaussian(m=1.0, V=1.0)), Message(Gaussian(m=2.0, V=-2.0)), nothing],
                                    Gaussian(m=3.0, V=-1.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=2.0, V=-2.0)), Message(Gaussian(m=3.0, V=1.0))],
                                    Gaussian(m=1.0, V=-1.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(Gaussian(m=1.0, V=-2.0)), nothing, Message(Gaussian(m=3.0, V=1.0))],
                                    Gaussian(m=2.0, V=-1.0))
            addition_node = AdditionNode()
            @fact_throws sumProductRule!(addition_node, 3, [Message(Gaussian(m=1.0, V=-1.0)), Message(Gaussian(m=2.0, V=-2.0)), nothing])
            @fact_throws sumProductRule!(addition_node, 1, [nothing, Message(Gaussian(m=1.0, V=-1.0)), Message(Gaussian(m=2.0, V=-2.0))])
        end
    end

    context("AdditionNode should provide sumProductRule! for combinations of Gaussian and Delta") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(Gaussian(m=1.0, V=0.5)), Message(Delta(3.0)), nothing],
                                Gaussian(m=4.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(Delta(3.0)), Message(Gaussian(m=1.0, V=0.5)), nothing],
                                Gaussian(m=4.0, V=0.5+tiny))
        #Backward message towards in1
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(Delta(3.0)), Message(Gaussian(m=1.0, V=0.5))],
                                Gaussian(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(Gaussian(m=1.0, V=0.5)), Message(Delta(3.0))],
                                Gaussian(m=2.0, V=0.5+tiny))
        # Backward message towards in2
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(Delta(3.0)), nothing, Message(Gaussian(m=1.0, V=0.5))],
                                Gaussian(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(Gaussian(m=1.0, V=0.5)), nothing, Message(Delta(3.0))],
                                Gaussian(m=2.0, V=0.5+tiny))
    end

    context("AdditionNode should provide sumProductRule! for MvGaussian") do
        context("MvGaussian with (m,V) parametrization") do
            mean = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussian(m=mean, V=variance)), Message(MvGaussian(m=mean, V=variance)), nothing],
                                    MvGaussian(m=[2.0, 4.0, 6.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussian(m=mean, V=variance)), Message(MvGaussian(m=mean, V=variance))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussian(m=mean, V=variance)), nothing, Message(MvGaussian(m=mean, V=variance))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("MvGaussian with (m,W) parametrization") do
            mean = collect(1.0:3.0)
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussian(m=mean, W=precision)), Message(MvGaussian(m=mean, W=precision)), nothing],
                                    MvGaussian(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussian(m=mean, W=precision)), Message(MvGaussian(m=mean, W=precision))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussian(m=mean, W=precision)), nothing, Message(MvGaussian(m=mean, W=precision))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end

        context("MvGaussian with (xi,V) parametrization") do
            xi = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussian(xi=xi, V=variance)), Message(MvGaussian(xi=xi, V=variance)), nothing],
                                    MvGaussian(xi=[1.0, 2.0, 3.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussian(xi=xi, V=variance)), Message(MvGaussian(xi=xi, V=variance))],
                                    MvGaussian(xi=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussian(xi=xi, V=variance)), nothing, Message(MvGaussian(xi=xi, V=variance))],
                                    MvGaussian(xi=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("MvGaussian with different parametrizations") do
            mean = collect(1.0:3.0)
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussian(m=mean, W=precision)), Message(MvGaussian(xi=precision*mean, V=inv(precision))), nothing],
                                    MvGaussian(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussian(m=mean, W=precision)), Message(MvGaussian(xi=precision*mean, V=inv(precision)))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussian(m=mean, W=precision)), nothing, Message(MvGaussian(xi=precision*mean, V=inv(precision)))],
                                    MvGaussian(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end
    end
end
