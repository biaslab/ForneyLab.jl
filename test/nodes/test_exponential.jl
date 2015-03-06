#####################
# Unit tests
#####################

facts("ExponentialNode unit tests") do
    context("ExponentialNode should initialize a ExponentialNode with 2 interfaces") do
        node = ExponentialNode()
        @fact typeof(node) => ExponentialNode
        @fact length(node.interfaces) => 2
        @fact node.in1 => node.interfaces[1]
        @fact node.out => node.interfaces[2]
    end

    context("ExponentialNode should pass messages") do
        # Forward message
        validateOutboundMessage(ExponentialNode(), 
                                2, 
                                [Message(GaussianDistribution(m=1.0, W=1.0)), nothing],
                                GammaDistribution(a=2.0, b=1.0/e))
        # Backward message
        validateOutboundMessage(ExponentialNode(), 
                                1, 
                                [nothing, Message(GammaDistribution(a=2.0, b=e))],
                                GaussianDistribution(m=-1.0, W=1.0))
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
