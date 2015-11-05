#####################
# Unit tests
#####################

facts("ExponentialNode unit tests") do
    context("ExponentialNode should initialize a ExponentialNode with 2 interfaces") do
        FactorGraph()
        ExponentialNode(id=:node)
        @fact typeof(n(:node)) --> ExponentialNode
        @fact length(n(:node).interfaces) --> 2
        @fact n(:node).i[:in] --> n(:node).interfaces[1]
        @fact n(:node).i[:out] --> n(:node).interfaces[2]
    end

    context("ExponentialNode should pass messages") do
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
end
