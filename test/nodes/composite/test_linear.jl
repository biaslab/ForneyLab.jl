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

    context("LinearCompositeNode should propagate a forward message to out") do
        validateOutboundMessage(LinearCompositeNode(form="precision"), 
                                5, 
                                GaussianDistribution, 
                                [Message(GaussianDistribution(m=1.0, W=1.0)), Message(0.5), Message(1.5), Message(1.0), nothing],
                                GaussianDistribution(m=2.0, W=0.8))
    end

    context("LinearCompositeNode should propagate a backward message to in1") do
        validateOutboundMessage(LinearCompositeNode(form="precision"), 
                                1, 
                                GaussianDistribution, 
                                [nothing, Message(0.5), Message(1.5), Message(1.0), Message(GaussianDistribution(m=2.0, W=0.8))],
                                GaussianDistribution(m=1.0, W=1/9))
    end

    context("LinearCompositeNode should propagate a backward variational message to in1") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 1, GaussianDistribution, nothing, GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [1.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to slope") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 2, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), nothing, GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [2.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to offset") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 3, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), nothing, InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=2.5, V=0.0))
        @fact msg.payload.m => [0.5]
    end

    context("LinearCompositeNode should propagate a backward variational message to noise") do
        validateOutboundMessage(LinearCompositeNode(), 
                                4, 
                                InverseGammaDistribution, 
                                [GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), nothing, GaussianDistribution(m=2.5, V=0.0)],
                                InverseGammaDistribution(a=-0.5, b=0.0))
    end

    context("LinearCompositeNode should propagate a forward variational message to out") do
        msg = ForneyLab.updateNodeMessage!(LinearCompositeNode(), 5, GaussianDistribution, GaussianDistribution(m=1.0, V=0.0), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), InverseGammaDistribution(a=10000.0, b=19998.0), nothing)
        @fact msg.payload.m => [2.5]
    end

    context("LinearCompositeNode should propagate a backward variational message to in1 under structured factorization") do
        validateOutboundMessage(LinearCompositeNode(), 
                                1, 
                                GaussianDistribution, 
                                [nothing, GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), GammaDistribution(a=0.1, b=0.1), Message(GaussianDistribution(m=2.5, V=0.0))],
                                GaussianDistribution(m=0.25, W=4.0))
    end

    context("LinearCompositeNode should propagate a backward variational message to slope under structured factorization") do
        dist_xy = GaussianDistribution(m=[1.0, 2.0], V=[1.0 0.5; 0.5 0.9])
        validateOutboundMessage(LinearCompositeNode(), 
                                2, 
                                GaussianDistribution, 
                                [dist_xy, nothing, GaussianDistribution(m=0.5, V=0.0), GammaDistribution(a=0.1, b=0.1), dist_xy],
                                GaussianDistribution(m=1.0, W=2.0))
    end

    context("LinearCompositeNode should propagate a backward variational message to offset under structured factorization") do
        dist_xy = GaussianDistribution(m=[1.0, 2.0], V=[1.0 0.5; 0.5 0.9])
        validateOutboundMessage(LinearCompositeNode(), 
                                3, 
                                GaussianDistribution, 
                                [dist_xy, GaussianDistribution(m=2.0, V=0.0), nothing, GammaDistribution(a=0.1, b=0.1), dist_xy],
                                GaussianDistribution(m=0.0, W=1.0))
    end

    context("LinearCompositeNode should propagate a backward variational message to noise under structured factorization") do
        dist_xy = GaussianDistribution(m=[1.0, 2.0], V=[1.0 0.5; 0.5 0.9])
        validateOutboundMessage(LinearCompositeNode(), 
                                4, 
                                GammaDistribution, 
                                [dist_xy, GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), nothing, dist_xy],
                                GammaDistribution(a=1.5, b=3.425))
    end

    context("LinearCompositeNode should propagate a forward variational message to out under structured factorization") do
        validateOutboundMessage(LinearCompositeNode(), 
                                5,
                                GaussianDistribution, 
                                [Message(GaussianDistribution(m=1.0, V=0.0)), GaussianDistribution(m=2.0, V=0.0), GaussianDistribution(m=0.5, V=0.0), GammaDistribution(a=0.1, b=0.1), nothing],
                                GaussianDistribution(m=0.5, W=1.0))
    end
end