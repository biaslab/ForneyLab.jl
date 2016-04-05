#####################
# Unit tests
#####################

facts("ExponentialNode unit tests") do
    context("ExponentialNode should initialize a ExponentialNode with 2 interfaces") do
        FactorGraph()
        ExponentialNode(id=:node)
        @fact length(n(:node).interfaces) --> 2
        @fact n(:node).i[:in] --> n(:node).interfaces[1]
        @fact n(:node).i[:out] --> n(:node).interfaces[2]
    end

    context("ExponentialNode should pass Gaussian messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(GaussianDistribution(m=2.0, V=3.0)), nothing],
                                LogNormalDistribution(m=2.0, s=3.0))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(LogNormalDistribution(m=2.0, s=3.0))],
                                GaussianDistribution(m=2.0, V=3.0))
    end

    context("ExponentialNode should provide approximate rules for Gaussian messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(GaussianDistribution(m=2.0, V=3.0)), nothing],
                                GammaDistribution(a=4.0/3.0, b=2.0/3.0),
                                sumProductRule!,
                                MomentMatching)
        a = ForneyLab.trigammaInverse(3.0)
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(GaussianDistribution(m=2.0, V=3.0)), nothing],
                                GammaDistribution(a=a, b=1/exp(2.0-digamma(a))),
                                sumProductRule!,
                                LogMomentMatching)
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(GammaDistribution(a=2.0, b=3.0))],
                                GaussianDistribution(m=digamma(2.0)-log(3.0), V=trigamma(2.0)),
                                sumProductRule!,
                                MomentMatching)
    end

    context("ExponentialNode should pass delta messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(DeltaDistribution(2.0)), nothing],
                                DeltaDistribution(exp(2.0)))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(DeltaDistribution(2.0))],
                                DeltaDistribution(log(2.0)))
    end

    context("ExponentialNode should pass multivariate Gaussian messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(MvGaussianDistribution(m=[1.0, 2.0], V=3.0*eye(2))), nothing],
                                MvLogNormalDistribution(m=[1.0, 2.0], S=3.0*eye(2)))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(MvLogNormalDistribution(m=[1.0, 2.0], S=3.0*eye(2)))],
                                MvGaussianDistribution(m=[1.0, 2.0], V=3.0*eye(2)))
    end

    context("ExponentialNode should pass multivariate delta messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(MvDeltaDistribution([1.0, 2.0])), nothing],
                                MvDeltaDistribution(exp([1.0, 2.0])))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(MvDeltaDistribution([1.0, 2.0]))],
                                MvDeltaDistribution(log([1.0, 2.0])))
    end
end
