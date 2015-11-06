#####################
# Unit tests
#####################

facts("EqualityNode unit tests") do
   context("EqualityNode() should initialize an EqualityNode with 3 interfaces") do
        FactorGraph()
        EqualityNode(id=:node)
        @fact typeof(n(:node)) --> EqualityNode
        @fact length(n(:node).interfaces) --> 3
        @fact ForneyLab.firstFreeInterface(n(:node)) --> n(:node).interfaces[1]
        for i=1:3
            Edge(n(:node), TerminalNode())
        end
        @fact_throws ForneyLab.firstFreeInterface(n(:node))
    end

    context("EqualityNode should provide sumProduct! for arbitrary message types") do
        # Equality constraint node should work for arbitraty messages, although not really useful.
        # Outbound message is equal to the inbound messages if not all inbound messages are equal.
        # Otherwise, the outbound message is Message(DeltaDistribution(0.0))

        # Equal scalars
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(DeltaDistribution(1.0)), Message(DeltaDistribution(1.0)), nothing],
                                DeltaDistribution(1.0))
        # Unequal scalars
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(DeltaDistribution(1.0)), Message(DeltaDistribution(1.1)), nothing],
                                DeltaDistribution(0.0))
    end

    context("EqualityNode should provide sumProduct! for GaussianDistribution") do
        inbound_dist = GaussianDistribution(xi=0.6, W=0.2)
        xi = inbound_dist.xi
        W = inbound_dist.W
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(xi=(xi + xi), W=(W + W)))
    end

    context("EqualityNode should provide sumProduct! for MvGaussianDistribution") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of sumProduct!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("MvGaussianDistribution with (m,V) parametrization") do
            mean = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = MvGaussianDistribution(m=mean, V=variance)
            W = inv(inbound_dist.V)
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("MvGaussianDistribution with (m,W) parametrization") do
            mean = collect(1.0:3.0)
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = MvGaussianDistribution(m=mean, W=precision)
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), W=(W+W)))
        end
        context("MvGaussianDistribution with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = MvGaussianDistribution(xi=inv(variance)*collect(1.0:3.0), V=variance)
            xi = inbound_dist.xi
            V = inbound_dist.V
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussianDistribution(xi=(xi + xi), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("MvGaussianDistribution with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = MvGaussianDistribution(xi=precision*collect(1.0:3.0), W=precision)
            xi = inbound_dist.xi
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussianDistribution(xi=(xi + xi), W=(W+W)))
        end
    end

    context("EqualityNode should provide sumProduct! for GammaDistribution") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(GammaDistribution()), Message(GammaDistribution()), nothing],
                                GammaDistribution(a=1.0, b=2.0))
    end

    context("EqualityNode should provide sumProduct! for InverseGammaDistribution") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(InverseGammaDistribution()), Message(InverseGammaDistribution()), nothing],
                                InverseGammaDistribution(a=7.0, b=4.0))
    end

    context("EqualityNode should provide sumProduct! for BetaDistribution") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(BetaDistribution(a=1.0, b=2.0)), Message(BetaDistribution(a=3.0, b=4.0)), nothing],
                                BetaDistribution(a=3.0, b=5.0))
    end

    context("EqualityNode should provide sumProduct! for combination of student's t and Gaussian distribution") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0)), Message(GaussianDistribution(m=0.0, V=1.0)), nothing],
                                GaussianDistribution(m=0.5, W=2.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(GaussianDistribution(m=0.0, V=1.0)), Message(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0)), nothing],
                                GaussianDistribution(m=0.5, W=2.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0)), nothing, Message(GaussianDistribution(m=0.0, V=1.0))],
                                GaussianDistribution(m=0.5, W=2.0))
    end

    context("EqualityNode should provide sumProduct! for combination of DeltaDistribution and a GaussianDistribution") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(DeltaDistribution(5.0)), Message(GaussianDistribution()), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(GaussianDistribution()), Message(DeltaDistribution(5.0)), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(DeltaDistribution(5.0)), nothing, Message(GaussianDistribution())],
                                DeltaDistribution(5.0))
    end

    context("EqualityNode should provide sumProduct! for combination of DeltaDistribution and a GammaDistribution") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(DeltaDistribution(5.0)), Message(GammaDistribution()), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(GammaDistribution()), Message(DeltaDistribution(5.0)), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(DeltaDistribution(5.0)), nothing, Message(GammaDistribution())],
                                DeltaDistribution(5.0))
    end

    context("EqualityNode should provide sumProduct! for combination of DeltaDistribution and a LogNormalDistribution") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(DeltaDistribution(5.0)), Message(LogNormalDistribution()), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(LogNormalDistribution()), Message(DeltaDistribution(5.0)), nothing],
                                DeltaDistribution(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(DeltaDistribution(5.0)), nothing, Message(LogNormalDistribution())],
                                DeltaDistribution(5.0))
    end
end
