#####################
# Unit tests
#####################

facts("EqualityNode unit tests") do
   context("EqualityNode() should initialize an EqualityNode with 3 interfaces") do
        FactorGraph()
        EqualityNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact ForneyLab.firstFreeInterface(n(:node)) --> n(:node).interfaces[1]
        for i=1:3
            Edge(n(:node), TerminalNode())
        end
        @fact_throws ForneyLab.firstFreeInterface(n(:node))
    end

    context("EqualityNode should provide sumProductRule! for Delta") do
        # Equality constraint node should work for arbitraty messages, although not really useful.
        # Outbound message is equal to the inbound messages if not all inbound messages are equal.
        # Otherwise, the outbound message is Message(Delta(0.0))

        # Equal scalars
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Delta(1.0)), Message(Delta(1.0)), nothing],
                                Delta(1.0))
    end

    context("EqualityNode should provide sumProductRule! for Gaussian") do
        context("Proper input distributions") do
            inbound_dist = Gaussian(xi=0.6, W=0.2)
            validateOutboundMessage(EqualityNode(),
                                    3,
                                    [Message(inbound_dist), Message(inbound_dist), nothing],
                                    Gaussian(xi=(0.6 + 0.6), W=(0.2 + 0.2)))
        end
        context("Improper input distributions") do
            inbound_dist = Gaussian(xi=0.6, W=0.2)
            validateOutboundMessage(EqualityNode(),
                                    3,
                                    [Message(Gaussian(xi=0.6, W=1.0)), Message(Gaussian(xi=0.6, W=-0.2)), nothing],
                                    Gaussian(xi=1.2, W=0.8))
        end
    end

    context("EqualityNode should provide sumProductRule! for MvGaussian") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of sumProductRule!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("MvGaussian with (m,V) parametrization") do
            mean = collect(1.0:3.0)
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = MvGaussian(m=mean, V=variance)
            W = inv(inbound_dist.V)
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussian(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("MvGaussian with (m,W) parametrization") do
            mean = collect(1.0:3.0)
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = MvGaussian(m=mean, W=precision)
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussian(m=(pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m)), W=(W+W)))
        end
        context("MvGaussian with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = MvGaussian(xi=inv(variance)*collect(1.0:3.0), V=variance)
            xi = inbound_dist.xi
            V = inbound_dist.V
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussian(xi=(xi + xi), V=(inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)))
        end
        context("MvGaussian with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = MvGaussian(xi=precision*collect(1.0:3.0), W=precision)
            xi = inbound_dist.xi
            W = inbound_dist.W
            validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(inbound_dist), Message(inbound_dist), nothing],
                                MvGaussian(xi=(xi + xi), W=(W+W)))
        end
    end

    context("EqualityNode should provide sumProductRule! for Gamma") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Gamma()), Message(Gamma()), nothing],
                                Gamma(a=1.0, b=2.0))
    end

    context("EqualityNode should provide sumProductRule! for InverseGamma") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(InverseGamma()), Message(InverseGamma()), nothing],
                                InverseGamma(a=7.0, b=4.0))
    end

    context("EqualityNode should provide sumProductRule! for Beta") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Beta(a=1.0, b=2.0)), Message(Beta(a=3.0, b=4.0)), nothing],
                                Beta(a=3.0, b=5.0))
    end

    context("EqualityNode should provide sumProductRule! for Bernoulli") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Bernoulli(0.2)), Message(Bernoulli(0.4)), nothing],
                                Bernoulli(1/7))
        validateOutboundMessage(EqualityNode(),
                                1,
                                [nothing, Message(Bernoulli(0.2)), Message(Bernoulli(0.4))],
                                Bernoulli(1/7))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(Bernoulli(0.2)), nothing, Message(Bernoulli(0.4))],
                                Bernoulli(1/7))
    end

    context("EqualityNode should provide sumProductRule! for Categorical") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Categorical([0.4;0.6])), Message(Categorical([0.4;0.6])), nothing],
                                Categorical([0.3076923076923077;0.6923076923076922]))
    end

    context("EqualityNode should provide approximate sumProductRule! for combination of student's t and Gaussian distribution") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(StudentsT(m=1.0, lambda=2.0, nu=4.0)), Message(Gaussian(m=0.0, V=1.0)), nothing, MomentMatching],
                                Gaussian(m=0.5, W=2.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Gaussian(m=0.0, V=1.0)), Message(StudentsT(m=1.0, lambda=2.0, nu=4.0)), nothing, MomentMatching],
                                Gaussian(m=0.5, W=2.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(StudentsT(m=1.0, lambda=2.0, nu=4.0)), nothing, Message(Gaussian(m=0.0, V=1.0)), MomentMatching],
                                Gaussian(m=0.5, W=2.0))
    end

    context("EqualityNode should provide sumProductRule! for combination of Delta and a Gaussian") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Delta(5.0)), Message(Gaussian()), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Gaussian()), Message(Delta(5.0)), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(Delta(5.0)), nothing, Message(Gaussian())],
                                Delta(5.0))
    end

    context("EqualityNode should provide sumProductRule! for combination of MvDelta and a MvGaussian") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(MvDelta([4.0, 5.0])), Message(MvGaussian(m=[1.0, 2.0], V=[1.0 0.5; 0.5 1.0])), nothing],
                                MvDelta([4.0, 5.0]))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(MvGaussian(m=[1.0, 2.0], V=[1.0 0.5; 0.5 1.0])), Message(MvDelta([4.0, 5.0])), nothing],
                                MvDelta([4.0, 5.0]))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(MvGaussian(m=[1.0, 2.0], V=[1.0 0.5; 0.5 1.0])), nothing, Message(MvDelta([4.0, 5.0]))],
                                MvDelta([4.0, 5.0]))
    end

    context("EqualityNode should provide sumProductRule! for combination of MatrixDelta and a Wishart") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(MatrixDelta([1.0 0.5; 0.5 1.0])), Message(Wishart(V=[2.0 1.0; 1.0 2.0], nu=3.0)), nothing],
                                MatrixDelta([1.0 0.5; 0.5 1.0]))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Wishart(V=[2.0 1.0; 1.0 2.0], nu=3.0)), Message(MatrixDelta([1.0 0.5; 0.5 1.0])), nothing],
                                MatrixDelta([1.0 0.5; 0.5 1.0]))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(Wishart(V=[2.0 1.0; 1.0 2.0], nu=3.0)), nothing, Message(MatrixDelta([1.0 0.5; 0.5 1.0]))],
                                MatrixDelta([1.0 0.5; 0.5 1.0]))
    end

    context("EqualityNode should provide sumProductRule! for combination of Delta and a Gamma") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Delta(5.0)), Message(Gamma()), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Gamma()), Message(Delta(5.0)), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(Delta(5.0)), nothing, Message(Gamma())],
                                Delta(5.0))
    end

    context("EqualityNode should provide sumProductRule! for combination of Delta and a LogNormal") do
        # Just test the original and a permutation of the arguments
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Delta(5.0)), Message(LogNormal()), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(LogNormal()), Message(Delta(5.0)), nothing],
                                Delta(5.0))
        validateOutboundMessage(EqualityNode(),
                                2,
                                [Message(Delta(5.0)), nothing, Message(LogNormal())],
                                Delta(5.0))
    end
end
