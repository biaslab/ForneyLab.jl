#####################
# Unit tests
#####################

facts("EqualityNode unit tests") do
   context("EqualityNode() should initialize an EqualityNode with 3 interfaces") do
        FactorGraph()
        node = EqualityNode()
        @fact typeof(node) => EqualityNode
        @fact length(node.interfaces) => 3
        @fact ForneyLab.firstFreeInterface(node) => node.interfaces[1]
        for i=1:3
            Edge(node, TerminalNode())
        end
        @fact_throws ForneyLab.firstFreeInterface(node)
    end

    context("EqualityNode should propagate an arbitrary message") do
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
        # Equal matrices
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(DeltaDistribution(ones(2,2))), Message(DeltaDistribution(ones(2,2))), nothing],
                                DeltaDistribution(ones(2,2)))
        # Unequal matrices (different values) should give zeros matrix
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(DeltaDistribution(ones(2,2))), Message(DeltaDistribution(4.0*ones(2,2))), nothing],
                                DeltaDistribution(zeros(2,2)))
    end

    context("EqualityNode should propagate a univariate GaussianDistribution") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of sumProduct!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("EqualityNode should propagate a univariate GaussianDistribution with (m,V) parametrization") do
            inbound_dist = GaussianDistribution(m=3.0, V=5.0)
            W = inv(inbound_dist.V)
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("EqualityNode should propagate a univariate GaussianDistribution with (m,W) parametrization") do
            inbound_dist = GaussianDistribution(m=3.0, W=0.2)
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), W=(W + W)))
        end
        context("EqualityNode should propagate a univariate GaussianDistribution with (xi,V) parametrization") do
            inbound_dist = GaussianDistribution(xi=0.6, V=5.0)
            xi = inbound_dist.xi
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(xi=(xi + xi), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("EqualityNode should propagate a univariate GaussianDistribution with (xi,W) parametrization") do
            inbound_dist = GaussianDistribution(xi=0.6, W=0.2)
            xi = inbound_dist.xi
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(xi=(xi + xi), W=(W + W)))
        end
    end

    context("EqualityNode should propagate a multivariate GaussianDistribution") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of sumProduct!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("EqualityNode should propagate a multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = GaussianDistribution(m=mean, V=variance)
            W = inv(inbound_dist.V)
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("EqualityNode should propagate a multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = GaussianDistribution(m=mean, W=precision)
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), W=(W+W)))
        end
        context("EqualityNode should propagate a multivariate GaussianDistribution with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = GaussianDistribution(xi=inv(variance)*[1.0:3.0], V=variance)
            xi = inbound_dist.xi
            V = inbound_dist.V
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(xi=(xi + xi), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("EqualityNode should propagate a multivariate GaussianDistribution with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = GaussianDistribution(xi=precision*[1.0:3.0], W=precision)
            xi = inbound_dist.xi
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                GaussianDistribution(xi=(xi + xi), W=(W+W)))
        end
    end

    context("EqualityNode should propagate a GammaMessage") do
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(GammaDistribution()), Message(GammaDistribution()), nothing],
                                GammaDistribution(a=1.0, b=2.0))
    end

    context("EqualityNode should propagate an InverseGammaMessage") do
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(InverseGammaDistribution()), Message(InverseGammaDistribution()), nothing],
                                InverseGammaDistribution(a=7.0, b=4.0))
    end

    context("EqualityNode should propagate a BetaMessage") do
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(BetaDistribution(a=1.0, b=2.0)), Message(BetaDistribution(a=3.0, b=4.0)), nothing],
                                BetaDistribution(a=3.0, b=5.0))
    end

    context("EqualityNode should propagate combination of student's t and Gaussian distribution") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(StudentsTDistribution()), Message(GaussianDistribution()), nothing],
                                GaussianDistribution(m=0.0, W=2.0))
        validateOutboundMessage(EqualityNode(), 
                                3, 
                                [Message(GaussianDistribution()), Message(StudentsTDistribution()), nothing],
                                GaussianDistribution(m=0.0, W=2.0))
        validateOutboundMessage(EqualityNode(), 
                                2, 
                                [Message(StudentsTDistribution()), nothing, Message(GaussianDistribution())],
                                GaussianDistribution(m=0.0, W=2.0))
    end

    context("EqualityNode should propagate combination of a delta and a Gaussian distribution") do
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

    context("EqualityNode should propagate combination of a delta and a gamma distribution") do
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
end
