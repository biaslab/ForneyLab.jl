#####################
# Unit tests
#####################

facts("GaussianNode unit tests") do
    context("Construction") do
        FactorGraph()

        testnode = GaussianNode(id=:node)
        @fact typeof(testnode) --> GaussianNode{Val{:mean},Val{:variance}}
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:mean] --> n(:node).interfaces[1]
        @fact n(:node).i[:variance] --> n(:node).interfaces[2]
        @fact n(:node).i[:out] --> n(:node).interfaces[3]

        testnode = GaussianNode(form=:precision)
        @fact typeof(testnode) --> GaussianNode{Val{:mean},Val{:precision}}
        @fact testnode.i[:precision] --> testnode.interfaces[2]

        testnode = GaussianNode(form=:variance, m=1.0)
        @fact typeof(testnode) --> GaussianNode{Val{:fixed_mean},Val{:variance}}
        @fact testnode.i[:variance] --> testnode.interfaces[1]

        testnode = GaussianNode(form=:variance, m=1.0, V=1.0)
        @fact typeof(testnode) --> GaussianNode{Val{:fixed_mean},Val{:fixed_variance}}
        @fact testnode.i[:out] --> testnode.interfaces[1]

        testnode = GaussianNode(form=:variance, V=1.0)
        @fact typeof(testnode) --> GaussianNode{Val{:mean},Val{:fixed_variance}}
        @fact testnode.i[:mean] --> testnode.interfaces[1]
        @fact testnode.i[:out] --> testnode.interfaces[2]

        testnode = GaussianNode(form=:log_variance)
        @fact typeof(testnode) --> GaussianNode{Val{:mean},Val{:log_variance}}
        @fact testnode.i[:log_variance] --> testnode.interfaces[2]
    end

    FactorGraph()

    context("Sum-product: no fixed parameters") do
        context("Sum-product message towards out") do
            validateOutboundMessage(GaussianNode(),
                                    3,
                                    [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("Sum-product message towards mean") do
            validateOutboundMessage(GaussianNode(),
                                    1,
                                    [nothing, Message(DeltaDistribution(0.5)), Message(DeltaDistribution(2.0))],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    1,
                                    [nothing, Message(DeltaDistribution(0.5)), Message(DeltaDistribution(2.0))],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("Sum-product message towards variance/precision") do
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

    context("Sum-product: fixed mean") do
        context("Sum-product message towards out") do
            validateOutboundMessage(GaussianNode(m=2.0),
                                    2,
                                    [Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(m=2.0; form=:precision),
                                    2,
                                    [Message(DeltaDistribution(0.5)), nothing],
                                    GaussianDistribution(m=2.0, W=0.5))
        end

        context("Sum-product message towards variance/precision") do
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

    context("Sum-product: fixed variance/precision") do
        context("Sum-product message towards out") do
            validateOutboundMessage(GaussianNode(V=0.5),
                                    2,
                                    [Message(DeltaDistribution(2.0)), nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
        end

        context("Sum-product message towards mean") do
            validateOutboundMessage(GaussianNode(V=0.5),
                                    1,
                                    [nothing, Message(DeltaDistribution(1.0))],
                                    GaussianDistribution(m=1.0, V=0.5))
        end
    end

    context("Sum-product: all parameters fixed") do
            validateOutboundMessage(GaussianNode(m=2.0, V=0.5),
                                    1,
                                    [nothing],
                                    GaussianDistribution(m=2.0, V=0.5))
    end

    context("Variational estimation") do
        context("GaussianNode with fixed variance should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(V=2.0),
                                    2,
                                    [GaussianDistribution(m=3.0, V=1.0), nothing],
                                    GaussianDistribution(m=3.0, V=2.0),
                                    ForneyLab.variationalRule!)
        end

        context("Naive variational implementation (mean field)") do
            context("GaussianNode should propagate a backward message to the mean") do
                # Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, GammaDistribution(a=3.0, b=1.0), GaussianDistribution(m=2.0, V=0.1)],
                                        GaussianDistribution(m=2.0, W=3.0),
                                        ForneyLab.variationalRule!)
                # Mv Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, WishartDistribution(V=[1.0 0.0; 0.0 2.0], nu=2.0), MvGaussianDistribution(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0])],
                                        MvGaussianDistribution(m=[1.0, 2.0], W=[2.0 0.0; 0.0 4.0]),
                                        ForneyLab.variationalRule!)
                # Log-Variance
                validateOutboundMessage(GaussianNode(form=:log_variance),
                                        1,
                                        [nothing, GaussianDistribution(m=3.0, V=2.0), GaussianDistribution(m=2.0, V=0.1)],
                                        GaussianDistribution(m=2.0, V=exp(4.0)),
                                        ForneyLab.variationalRule!)
                # Moment
                validateOutboundMessage(GaussianNode(),
                                        1,
                                        [nothing, InverseGammaDistribution(a=3.0, b=1.0), GaussianDistribution(m=2.0, V=0.1)],
                                        GaussianDistribution(m=2.0, V=4.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a backward message to the variance, precision or log-precision") do
                # Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [GaussianDistribution(m=4.0, W=2.0), nothing, GaussianDistribution(m=2.0, V=0.1)],
                                        GammaDistribution(a=1.5, b=2.3),
                                        ForneyLab.variationalRule!)
                # MvPrecision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [MvGaussianDistribution(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0]), nothing, MvGaussianDistribution(m=[3.0, 4.0], V=[1.0 0.0; 0.0 2.0])],
                                        WishartDistribution(V=[0.25 -0.125; -0.125 0.1875], nu=4.0),
                                        ForneyLab.variationalRule!)
                # Log-Variance
                validateOutboundMessage(GaussianNode(form=:log_variance),
                                        2,
                                        [GaussianDistribution(m=3.0, V=2.0), nothing, GaussianDistribution(m=2.0, V=1.0)],
                                        GaussianDistribution(m=log(4.0), V=2.0),
                                        ForneyLab.variationalRule!)
                # Moment
                validateOutboundMessage(GaussianNode(form=:variance),
                                        2,
                                        [GaussianDistribution(m=4.0, V=1.0), nothing, GaussianDistribution(m=2.0, V=0.1)],
                                        InverseGammaDistribution(a=-0.5, b=2.55),
                                        ForneyLab.variationalRule!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        1,
                                        [nothing, GaussianDistribution(m=2.0, V=0.5)],
                                        GammaDistribution(a=1.5, b=0.75),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=1.0; form=:log_variance),
                                        1,
                                        [nothing, GaussianDistribution(m=3.0, V=2.0)],
                                        GaussianDistribution(m=log(6.0), V=2.0),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=4.0; form=:variance),
                                        1,
                                        [nothing, GaussianDistribution(m=2.0, V=0.1)],
                                        InverseGammaDistribution(a=-0.5, b=2.05),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a forward message") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [GaussianDistribution(), GammaDistribution(), nothing],
                                        GaussianDistribution(m=0.0, W=1.0),
                                        ForneyLab.variationalRule!)
                # Mv Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [MvGaussianDistribution(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0]), WishartDistribution(V=[1.0 0.0; 0.0 2.0], nu=2.0), nothing],
                                        MvGaussianDistribution(m=[1.0, 2.0], W=[2.0 0.0; 0.0 4.0]),
                                        ForneyLab.variationalRule!)
                # Log variance
                validateOutboundMessage(GaussianNode(form=:log_variance),
                                        3,
                                        [GaussianDistribution(m=2.0, V=0.1), GaussianDistribution(m=3.0, V=2.0), nothing],
                                        GaussianDistribution(m=2.0, V=exp(4.0)),
                                        ForneyLab.variationalRule!)
                # Variance
                validateOutboundMessage(GaussianNode(form=:variance),
                                        3,
                                        [GaussianDistribution(), InverseGammaDistribution(a=2.0, b=1.0), nothing],
                                        GaussianDistribution(m=0.0, V=1.0),
                                        ForneyLab.variationalRule!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        2,
                                        [GammaDistribution(a=1.5, b=0.5), nothing],
                                        GaussianDistribution(m=1.0, W=3.0),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=1.0; form=:log_variance),
                                        2,
                                        [GaussianDistribution(m=3.0, V=2.0), nothing],
                                        GaussianDistribution(m=1.0, V=exp(4.0)),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=1.0; form=:variance),
                                        2,
                                        [InverseGammaDistribution(a=2.0, b=1.0), nothing],
                                        GaussianDistribution(m=1.0, V=1.0),
                                        ForneyLab.variationalRule!)
            end
        end

        context("Structured variational implementation") do
            context("GaussianNode should propagate a backward message to the mean") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, Message(GammaDistribution(a=3.0, b=1.0)), GaussianDistribution(m=2.0, W=10.0)],
                                        StudentsTDistribution(m=2.0, lambda=60.0/21.0, nu=6.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a backward message to the precision") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [Message(GaussianDistribution(m=4.0, W=2.0)), nothing, GaussianDistribution(m=2.0, W=10.0)],
                                        GammaDistribution(a=1.5, b=41.0/20.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a forward message") do
                node_marg = NormalGammaDistribution(m=4.0, beta=2.0, a=3.0, b=1.0)
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [node_marg, node_marg, nothing],
                                        GaussianDistribution(m=4.0, W=3.0),
                                        ForneyLab.variationalRule!)
            end
        end
    end
end
