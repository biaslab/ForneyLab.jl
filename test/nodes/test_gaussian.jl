#####################
# Unit tests
#####################

facts("GaussianNode unit tests") do
    context("GaussianNode() should initialize a GaussianNode with 3 interfaces") do
        FactorGraph()
        GaussianNode(id=:node)
        @fact typeof(n(:node)) => GaussianNode
        @fact length(n(:node).interfaces) => 3
        @fact n(:node).i[:mean] => n(:node).interfaces[1]
        @fact n(:node).i[:variance] => n(:node).interfaces[2]
        @fact n(:node).i[:out] => n(:node).interfaces[3]
    end

    context("GaussianNode() should initialize a GaussianNode with precision parametrization") do
        FactorGraph()
        GaussianNode(form=:precision, id=:node)
        @fact n(:node).i[:mean] => n(:node).interfaces[1]
        @fact n(:node).i[:precision] => n(:node).interfaces[2]
        @fact n(:node).i[:out] => n(:node).interfaces[3]
    end

    FactorGraph()

    context("GaussianNode() should handle fixed mean") do
        context("GaussianNode with fixed mean should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(m=2.0),
                                    2,
                                    [Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(m=2.0; form=:precision),
                                    2,
                                    [Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("GaussianNode with fixed mean should propagate a backward message to the variance") do
            validateOutboundMessage(GaussianNode(m=2.0),
                                    1,
                                    [nothing, Message(DeltaDistribution(1.0))],
                                    InverseGammaDistribution(a=-0.5, b=0.5))
            validateOutboundMessage(GaussianNode(m=2.0; form=:precision),
                                    1,
                                    [nothing, Message(DeltaDistribution(1.0))],
                                    GammaDistribution(a=1.5, b=0.5))
        end
    end

    context("GaussianNode() should handle fixed mean and variance") do
            validateOutboundMessage(GaussianNode(m=2.0, V=0.5),
                                    1,
                                    [nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
    end

    context("GaussianNode() should handle fixed variance for mean field") do
        context("GaussianNode with fixed variance should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(V=2.0),
                                    2,
                                    [GaussianDistribution(m=3.0, V=1.0), nothing],
                                    GaussianDistribution(m=3.0, V=2.0),
                                    ForneyLab.vmp!)
        end

        context("GaussianNode with fixed variance should propagate a backward message to the variance") do
            validateOutboundMessage(GaussianNode(V=2.0),
                                    1,
                                    [nothing, GaussianDistribution(m=3.0, V=1.0)],
                                    GaussianDistribution(m=3.0, V=2.0),
                                    ForneyLab.vmp!)
        end
    end

    context("Point estimates of y and m, so no approximation is required.") do
        context("GaussianNode should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(),
                                    3,
                                    [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("GaussianNode should propagate a backward message to the mean") do
            validateOutboundMessage(GaussianNode(),
                                    1,
                                    [nothing, Message(DeltaDistribution(0.5)), Message(DeltaDistribution(2.0))],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    1,
                                    [nothing, Message(DeltaDistribution(0.5)), Message(DeltaDistribution(2.0))],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("GaussianNode should propagate a backward message to the variance/precision") do
            validateOutboundMessage(GaussianNode(),
                                    2,
                                    [Message(DeltaDistribution(2.0)), nothing, Message(DeltaDistribution(1.0))],
                                    InverseGammaDistribution(a=-0.5, b=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    2,
                                    [Message(DeltaDistribution(2.0)), nothing, Message(DeltaDistribution(1.0))],
                                    GammaDistribution(a=1.5, b=0.5))
        end
    end

    context("Variational estimation") do
        context("Naive variational implementation (mean field)") do
            context("GaussianNode should propagate a backward message to the mean") do
                # Standard
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, GammaDistribution(a=3.0, b=1.0), GaussianDistribution(m=2.0, V=0.1)],
                                        GaussianDistribution(m=2.0, W=3.0),
                                        ForneyLab.vmp!)
                # Inverse
                validateOutboundMessage(GaussianNode(),
                                        1,
                                        [nothing, InverseGammaDistribution(a=3.0, b=1.0), GaussianDistribution(m=2.0, V=0.1)],
                                        GaussianDistribution(m=2.0, V=4.0),
                                        ForneyLab.vmp!)
            end

            context("GaussianNode should propagate a backward message to the variance or precision") do
                # Standard
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [GaussianDistribution(m=4.0, W=2.0), nothing, GaussianDistribution(m=2.0, V=0.1)],
                                        GammaDistribution(a=1.5, b=2.3),
                                        ForneyLab.vmp!)
                # Inverse
                validateOutboundMessage(GaussianNode(),
                                        2,
                                        [GaussianDistribution(m=4.0, V=1.0), nothing, GaussianDistribution(m=2.0, V=0.1)],
                                        InverseGammaDistribution(a=-0.5, b=2.55),
                                        ForneyLab.vmp!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        1,
                                        [nothing, GaussianDistribution(m=2.0, V=0.5)],
                                        GammaDistribution(a=1.5, b=0.75),
                                        ForneyLab.vmp!)
            end

            context("GaussianNode should propagate a forward message") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [GaussianDistribution(), GammaDistribution(), nothing],
                                        GaussianDistribution(m=0.0, W=1.0),
                                        ForneyLab.vmp!)
                validateOutboundMessage(GaussianNode(form=:moment),
                                        3,
                                        [GaussianDistribution(), InverseGammaDistribution(a=2.0, b=1.0), nothing],
                                        GaussianDistribution(m=0.0, V=1.0),
                                        ForneyLab.vmp!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        2,
                                        [GammaDistribution(a=1.5, b=0.5), nothing],
                                        GaussianDistribution(m=1.0, W=3.0),
                                        ForneyLab.vmp!)
            end
        end

        context("Structured variational implementation") do
            context("GaussianNode should propagate a backward message to the mean") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, Message(GammaDistribution(a=3.0, b=1.0)), GaussianDistribution(m=2.0, W=10.0)],
                                        StudentsTDistribution(m=2.0, lambda=60.0/21.0, nu=6.0),
                                        ForneyLab.vmp!)
            end

            context("GaussianNode should propagate a backward message to the precision") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [Message(GaussianDistribution(m=4.0, W=2.0)), nothing, GaussianDistribution(m=2.0, W=10.0)],
                                        GammaDistribution(a=1.5, b=41.0/20.0),
                                        ForneyLab.vmp!)
            end

            context("GaussianNode should propagate a forward message") do
                node_marg = NormalGammaDistribution(m=4.0, beta=2.0, a=3.0, b=1.0)
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [node_marg, node_marg, nothing],
                                        GaussianDistribution(m=4.0, W=3.0),
                                        ForneyLab.vmp!)
            end
        end
    end
end