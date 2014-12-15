#####################
# Unit tests
#####################

facts("LogarithmicNode unit tests") do
    context("LogarithmicNode should initialize a LogarithmicNode with 2 interfaces") do
        node = LogarithmicNode()
        @fact typeof(node) => LogarithmicNode
        @fact length(node.interfaces) => 2
        @fact node.in1 => node.interfaces[1]
        @fact node.out => node.interfaces[2]
    end

    context("LogarithmicNode should pass messages") do
        # Forward message
        validateOutboundMessage(LogarithmicNode(), 
                                2, 
                                [Message(GaussianDistribution(m=1.0, W=1.0)), nothing],
                                GammaDistribution(a=2.0, b=1.0/e))
        # Backward message
        validateOutboundMessage(LogarithmicNode(), 
                                1, 
                                [nothing, Message(GammaDistribution(a=2.0, b=e))],
                                GaussianDistribution(m=-1.0, W=1.0))
    end
end
