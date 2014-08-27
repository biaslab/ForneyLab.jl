#####################
# Unit tests
#####################

facts("GaussianNode unit tests") do
    context("GaussianNode() should initialize a GaussianNode with 3 interfaces") do
        node = GaussianNode()
        @fact typeof(node) => GaussianNode
        @fact length(node.interfaces) => 3
        @fact node.mean => node.interfaces[1]
        @fact node.variance => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact node.variational => false # default variational to false
    end

    context("GaussianNode() should initialize a GaussianNode with precision parametrization") do
        node = GaussianNode(form="precision")
        @fact node.mean => node.interfaces[1]
        @fact node.precision => node.interfaces[2]
        @fact node.out => node.interfaces[3]
    end

    context("GaussianNode() should handle fixed parameters") do
        # Fix mean
        node = GaussianNode(m=GaussianDistribution())
        @fact typeof(node.mean.partner.node) => ForneyLab.ClampNode
        @fact node.mean.partner.message.payload => GaussianDistribution()
        # Fix variance
        node = GaussianNode(V=InverseGammaDistribution())
        @fact typeof(node.variance.partner.node) => ForneyLab.ClampNode
        @fact node.variance.partner.message.payload => InverseGammaDistribution()
        # Fix precision
        node = GaussianNode(form="precision", W=GammaDistribution())
        @fact typeof(node.precision.partner.node) => ForneyLab.ClampNode
        @fact node.precision.partner.message.payload => GammaDistribution()
        # Fix mean and variance
        node = GaussianNode(m=GaussianDistribution(), V=InverseGammaDistribution())
        @fact typeof(node.mean.partner.node) => ForneyLab.ClampNode
        @fact typeof(node.variance.partner.node) => ForneyLab.ClampNode
    end

    context("Point estimates of y and m, so no approximation is required.") do
        context("GaussianNode should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(2.0), Message(InverseGammaDistribution(a=3.0, b=1.0)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
        end

        context("GaussianNode should propagate a backward message to the mean") do
            validateOutboundMessage(GaussianNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(InverseGammaDistribution(a=3.0, b=1.0)), Message(2.0)],
                                    GaussianDistribution(m=2.0, V=0.5))
        end

        context("GaussianNode should propagate a backward message to the variance") do
            validateOutboundMessage(GaussianNode(), 
                                    2, 
                                    InverseGammaDistribution, 
                                    [Message(2.0), nothing, Message(1.0)],
                                    InverseGammaDistribution(a=-0.5, b=0.5))
        end
    end

    context("Variational estimation") do
        context("Naive variational implementation (mean field)") do
            context("GaussianNode should propagate a backward message to the mean") do
                # Standard
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        1, 
                                        GaussianDistribution, 
                                        [nothing, GammaDistribution(a=3.0, b=1.0), 2.0],
                                        GaussianDistribution(m=2.0, W=3.0))
                # Inverse
                validateOutboundMessage(GaussianNode(), 
                                        1, 
                                        GaussianDistribution, 
                                        [nothing, InverseGammaDistribution(a=3.0, b=1.0), 2.0],
                                        GaussianDistribution(m=2.0, V=4.0))
            end

            context("GaussianNode should propagate a backward message to the variance or precision") do
                # Standard
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        2, 
                                        GammaDistribution, 
                                        [GaussianDistribution(m=4.0, W=2.0), nothing, 2.0],
                                        GammaDistribution(a=1.5, b=2.25))
                # Inverse
                validateOutboundMessage(GaussianNode(), 
                                        2, 
                                        InverseGammaDistribution, 
                                        [GaussianDistribution(m=4.0, V=1.0), nothing, 2.0],
                                        InverseGammaDistribution(a=-0.5, b=2.5))
            end

            context("GaussianNode should propagate a forward message") do
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        3, 
                                        GaussianDistribution, 
                                        [GaussianDistribution(), GammaDistribution(), nothing],
                                        GaussianDistribution(m=0.0, W=1.0))
                validateOutboundMessage(GaussianNode(form="moment"), 
                                        3, 
                                        GaussianDistribution, 
                                        [GaussianDistribution(), InverseGammaDistribution(a=2.0, b=1.0), nothing],
                                        GaussianDistribution(m=0.0, V=1.0))
            end
        end

        context("Structured variational implementation") do
            context("GaussianNode should propagate a backward message to the mean") do
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        1, 
                                        StudentsTDistribution, 
                                        [nothing, Message(GammaDistribution(a=3.0, b=1.0)), GaussianDistribution(m=2.0, W=10.0)],
                                        StudentsTDistribution(m=2.0, W=60.0/21.0, nu=6.0))
            end

            context("GaussianNode should propagate a backward message to the precision") do
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        2, 
                                        GammaDistribution, 
                                        [Message(GaussianDistribution(m=4.0, W=2.0)), nothing, GaussianDistribution(m=2.0, W=10.0)],
                                        GammaDistribution(a=1.5, b=41.0/20.0))
            end

            context("GaussianNode should propagate a forward message") do
                validateOutboundMessage(GaussianNode(form="precision"), 
                                        3, 
                                        GaussianDistribution, 
                                        [NormalGammaDistribution(m=4.0, beta=2.0, a=3.0, b=1.0), nothing],
                                        GaussianDistribution(m=4.0, W=3.0))
            end
        end
    end

    context("Marginal calculation for the GaussianNode") do
        node = GaussianNode(form="precision")
        m_edge = Edge(MockNode(Message(GaussianDistribution())).out, node.mean)
        gam_edge = Edge(MockNode(Message(GammaDistribution())).out, node.precision, GammaDistribution)
        y_edge = Edge(node.out, MockNode(Message(GaussianDistribution(m=0.0, W=10.0))).out)
        y_edge.marginal = GaussianDistribution(m=0.0, W=10.0)
        calculateMarginal!(node)
        @fact node.marginal => NormalGammaDistribution(m=0.0, beta=10000.0, a=1.5, b=21.0/20.0) 
    end    
end