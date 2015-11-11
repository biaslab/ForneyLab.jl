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

    context("AdditionNode should provide sumProduct! for DeltaDistribution{Float64}") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0)), nothing],
                                DeltaDistribution(5.0))
        # Backward message
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))],
                                DeltaDistribution(1.0))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(DeltaDistribution(2.0)), nothing, Message(DeltaDistribution(3.0))],
                                DeltaDistribution(1.0))
    end

    context("AdditionNode should provide sumProduct! for MvDeltaDistribution{Float64}") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(MvDeltaDistribution([1.0, 2.0])), Message(MvDeltaDistribution([3.0, 4.0])), nothing],
                                MvDeltaDistribution([4.0, 6.0]))
        # Backward message
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(MvDeltaDistribution([1.0, 2.0])), Message(MvDeltaDistribution([3.0, 4.0]))],
                                MvDeltaDistribution([2.0, 2.0]))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(MvDeltaDistribution([1.0, 2.0])), nothing, Message(MvDeltaDistribution([3.0, 4.0]))],
                                MvDeltaDistribution([2.0, 2.0]))
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should provide sumProduct! for GaussianDistribution") do
        context("GaussianDistribution with (m,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0)), nothing],
                                    GaussianDistribution(m=4.0, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), nothing, Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))
        end

        context("GaussianDistribution with (m,W) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
        end

        context("GaussianDistribution with (xi,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0)), nothing],
                                    GaussianDistribution(xi=14/6, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), nothing, Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))
        end

        context("GaussianDistribution with different parametrizations") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
        end
        context("Support for improper GaussianDistributions") do
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(GaussianDistribution(m=1.0, V=1.0)), Message(GaussianDistribution(m=2.0, V=-2.0)), nothing],
                                    GaussianDistribution(m=3.0, V=-1.0))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(GaussianDistribution(m=2.0, V=-2.0)), Message(GaussianDistribution(m=3.0, V=1.0))],
                                    GaussianDistribution(m=1.0, V=-1.0))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(GaussianDistribution(m=1.0, V=-2.0)), nothing, Message(GaussianDistribution(m=3.0, V=1.0))],
                                    GaussianDistribution(m=2.0, V=-1.0))
            addition_node = AdditionNode()
            @fact_throws sumProduct!(addition_node, 3, [Message(GaussianDistribution(m=1.0, V=-1.0)), Message(GaussianDistribution(m=2.0, V=-2.0)), nothing])
            @fact_throws sumProduct!(addition_node, 1, [nothing, Message(GaussianDistribution(m=1.0, V=-1.0)), Message(GaussianDistribution(m=2.0, V=-2.0))])
        end
    end

    context("AdditionNode should provide sumProduct! for combinations of GaussianDistribution and DeltaDistribution") do
        # Forward message
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(GaussianDistribution(m=1.0, V=0.5)), Message(DeltaDistribution(3.0)), nothing],
                                GaussianDistribution(m=4.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                3,
                                [Message(DeltaDistribution(3.0)), Message(GaussianDistribution(m=1.0, V=0.5)), nothing],
                                GaussianDistribution(m=4.0, V=0.5+tiny))
        #Backward message towards in1
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(DeltaDistribution(3.0)), Message(GaussianDistribution(m=1.0, V=0.5))],
                                GaussianDistribution(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                1,
                                [nothing, Message(GaussianDistribution(m=1.0, V=0.5)), Message(DeltaDistribution(3.0))],
                                GaussianDistribution(m=2.0, V=0.5+tiny))
        # Backward message towards in2
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(DeltaDistribution(3.0)), nothing, Message(GaussianDistribution(m=1.0, V=0.5))],
                                GaussianDistribution(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(),
                                2,
                                [Message(GaussianDistribution(m=1.0, V=0.5)), nothing, Message(DeltaDistribution(3.0))],
                                GaussianDistribution(m=2.0, V=0.5+tiny))
    end

    context("AdditionNode should provide sumProduct! for MvGaussianDistribution") do
        context("MvGaussianDistribution with (m,V) parametrization") do
            mean = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussianDistribution(m=mean, V=variance)), Message(MvGaussianDistribution(m=mean, V=variance)), nothing],
                                    MvGaussianDistribution(m=[2.0, 4.0, 6.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussianDistribution(m=mean, V=variance)), Message(MvGaussianDistribution(m=mean, V=variance))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussianDistribution(m=mean, V=variance)), nothing, Message(MvGaussianDistribution(m=mean, V=variance))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("MvGaussianDistribution with (m,W) parametrization") do
            mean = collect(1.0:3.0)
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussianDistribution(m=mean, W=precision)), Message(MvGaussianDistribution(m=mean, W=precision)), nothing],
                                    MvGaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussianDistribution(m=mean, W=precision)), Message(MvGaussianDistribution(m=mean, W=precision))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussianDistribution(m=mean, W=precision)), nothing, Message(MvGaussianDistribution(m=mean, W=precision))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end

        context("MvGaussianDistribution with (xi,V) parametrization") do
            xi = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussianDistribution(xi=xi, V=variance)), Message(MvGaussianDistribution(xi=xi, V=variance)), nothing],
                                    MvGaussianDistribution(xi=[1.0, 2.0, 3.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussianDistribution(xi=xi, V=variance)), Message(MvGaussianDistribution(xi=xi, V=variance))],
                                    MvGaussianDistribution(xi=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussianDistribution(xi=xi, V=variance)), nothing, Message(MvGaussianDistribution(xi=xi, V=variance))],
                                    MvGaussianDistribution(xi=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("MvGaussianDistribution with different parametrizations") do
            mean = collect(1.0:3.0)
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(),
                                    3,
                                    [Message(MvGaussianDistribution(m=mean, W=precision)), Message(MvGaussianDistribution(xi=precision*mean, V=inv(precision))), nothing],
                                    MvGaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(),
                                    1,
                                    [nothing, Message(MvGaussianDistribution(m=mean, W=precision)), Message(MvGaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(),
                                    2,
                                    [Message(MvGaussianDistribution(m=mean, W=precision)), nothing, Message(MvGaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    MvGaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end
    end
end
