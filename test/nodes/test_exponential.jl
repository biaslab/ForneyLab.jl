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
                                [Message(Gaussian(m=2.0, V=3.0)), nothing],
                                LogNormal(m=2.0, s=3.0))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(LogNormal(m=2.0, s=3.0))],
                                Gaussian(m=2.0, V=3.0))
    end

    context("ExponentialNode should provide approximate rules for Gaussian messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(Gaussian(m=2.0, V=3.0)), nothing],
                                Gamma(a=4.0/3.0, b=2.0/3.0),
                                sumProductRule!,
                                MomentMatching)
        a = ForneyLab.trigammaInverse(3.0)
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(Gaussian(m=2.0, V=3.0)), nothing],
                                Gamma(a=a, b=1/exp(2.0-digamma(a))),
                                sumProductRule!,
                                LogMomentMatching)
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(Gamma(a=2.0, b=3.0))],
                                Gaussian(m=digamma(2.0)-log(3.0), V=trigamma(2.0)),
                                sumProductRule!,
                                MomentMatching)
    end

    context("ExponentialNode should pass delta messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(Delta(2.0)), nothing],
                                Delta(exp(2.0)))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(Delta(2.0))],
                                Delta(log(2.0)))
    end

    context("ExponentialNode should pass multivariate Gaussian messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(MvGaussian(m=[1.0, 2.0], V=3.0*eye(2))), nothing],
                                MvLogNormal(m=[1.0, 2.0], S=3.0*eye(2)))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(MvLogNormal(m=[1.0, 2.0], S=3.0*eye(2)))],
                                MvGaussian(m=[1.0, 2.0], V=3.0*eye(2)))
    end

    context("ExponentialNode should pass multivariate delta messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(),
                                2,
                                [Message(MvDelta([1.0, 2.0])), nothing],
                                MvDelta(exp([1.0, 2.0])))
        # Backward message
        validateOutboundMessage(ExponentialNode(),
                                1,
                                [nothing, Message(MvDelta([1.0, 2.0]))],
                                MvDelta(log([1.0, 2.0])))
    end
end
