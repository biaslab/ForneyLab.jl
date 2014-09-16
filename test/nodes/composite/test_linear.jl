#####################
# Unit tests
#####################

facts("LinearCompositeNode unit tests") do
    context("LinearCompositeNode() should initialize a LinearCompositeNode with 5 interfaces") do
        FactorGraph()
        node = LinearCompositeNode()
        @fact typeof(node) => LinearCompositeNode
        @fact length(node.interfaces) => 5
        @fact node.in1 => node.interfaces[1]
        @fact node.slope => node.interfaces[2]
        @fact node.offset => node.interfaces[3]
        @fact node.noise => node.interfaces[4]
        @fact node.out => node.interfaces[5]
        @fact node.use_composite_update_rules => true # default use_composite_update_rules to true
    end

    FactorGraph()

    context("LinearCompositeNode should propagate a backward variational message to in1") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 1, GaussianDistribution, uninformative(GaussianDistribution), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [1.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to slope") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 2, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), uninformative(GaussianDistribution), GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [2.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to offset") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 3, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), uninformative(GaussianDistribution), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [0.5]
    end

    context("LinearCompositeNode should propagate a backward variational message to noise") do
        validateOutboundMessage(LinearCompositeNode(), 
                                4, 
                                InverseGammaDistribution, 
                                [GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), uninformative(InverseGammaDistribution), GaussianDistribution(m=2.5, V=0.0)],
                                InverseGammaDistribution(a=-0.5, b=0.0))
    end

    context("LinearCompositeNode should propagate a forward variational message to out") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 5, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), uninformative(GaussianDistribution))
        @fact msg.payload.m => [2.5]
    end
end