#####################
# Unit tests
#####################

facts("SigmoidCompositeNode unit tests") do
    context("SigmoidCompositeNode() should initialize a SigmoidCompositeNode with 2 interfaces") do
        FactorGraph()
        node = SigmoidCompositeNode()
        @fact typeof(node) => SigmoidCompositeNode
        @fact length(node.interfaces) => 2
        @fact node.in1 => node.interfaces[1]
        @fact node.out => node.interfaces[2]
        @fact typeof(node.a) => Float64
        @fact typeof(node.b) => Float64
        @fact typeof(node.gamma) => Float64
    end

    context("SigmoidCompositeNode should propagate a forward sumproduct message to out") do
        validateOutboundMessage(SigmoidCompositeNode(a=0.5, b=1.5, gamma=2.0), 
                                2, 
                                [Message(GaussianDistribution(m=1.0, W=1.5)), nothing],
                                BetaDistribution(a=7.519299915529833, b=9.244388473653649))
    end

    context("SigmoidCompositeNode should propagate a backward sumproduct message to in1") do
        validateOutboundMessage(SigmoidCompositeNode(a=0.5, b=1.5, gamma=2.0), 
                                1, 
                                [nothing, Message(BetaDistribution(a=2.0, b=3.0))],
                                GaussianDistribution(m=-0.9241962407465937, W=0.19999999999999996))
    end

    context("SigmoidCompositeNode should propagate a forward variational message to out") do
        validateOutboundMessage(SigmoidCompositeNode(a=0.5, b=1.5, gamma=2.0), 
                                2, 
                                [GaussianDistribution(m=1.0, W=1.5), nothing],
                                BetaDistribution(a=11.455634179939883, b=36.32586116628848),
                                ForneyLab.vmp!)
    end

    context("SigmoidCompositeNode should propagate a backward variational message to in1") do
        validateOutboundMessage(SigmoidCompositeNode(a=0.5, b=1.5, gamma=2.0), 
                                1, 
                                [nothing, BetaDistribution(a=2.0, b=3.0)],
                                GaussianDistribution(m=-0.7324081924454064, W=0.2592),
                                ForneyLab.vmp!)
    end
end