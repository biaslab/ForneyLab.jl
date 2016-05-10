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
    end

    FactorGraph()

    context("Sum-product: no fixed parameters") do
        context("Sum-product message towards out") do
            # Univariate
            validateOutboundMessage(GaussianNode(),
                                    3,
                                    [Message(Delta(2.0)), Message(Delta(0.5)), nothing],
                                    Gaussian(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(Delta(2.0)), Message(Delta(0.5)), nothing],
                                    Gaussian(m=2.0, W=0.5))
            validateOutboundMessage(GaussianNode(),
                                    3,
                                    [Message(Gaussian(m=2.0, V=1.0)), Message(Delta(0.5)), nothing],
                                    Gaussian(m=2.0, V=1.5))

            # Multivariate
            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(MvDelta(2.0*ones(2))), Message(MatrixDelta(0.5*eye(2))), nothing],
                                    MvGaussian(m=2.0*ones(2), W=0.5*eye(2)))

            validateOutboundMessage(GaussianNode(),
                                    3,
                                    [Message(MvDelta(2.0*ones(2))), Message(MatrixDelta(0.5*diageye(2))), nothing],
                                    MvGaussian(m=2.0*ones(2), V=0.5*diageye(2)))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(MvDelta(2.0*ones(2))), Message(MatrixDelta(0.5*diageye(2))), nothing],
                                    MvGaussian(m=2.0*ones(2), W=0.5*diageye(2)))

            validateOutboundMessage(GaussianNode(form=:precision),
                                    3,
                                    [Message(MvGaussian(m=2.0*ones(2), V=1.0*diageye(2))), Message(MatrixDelta(0.5*eye(2))), nothing],
                                    MvGaussian(m=2.0*ones(2), V=3.0*eye(2)))
        end

        context("Sum-product message towards mean") do
            validateOutboundMessage(GaussianNode(),
                                    1,
                                    [nothing, Message(Delta(0.5)), Message(Delta(2.0))],
                                    Gaussian(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    1,
                                    [nothing, Message(Delta(0.5)), Message(Delta(2.0))],
                                    Gaussian(m=2.0, W=0.5))

            validateOutboundMessage(GaussianNode(),
                                    1,
                                    [nothing, Message(Delta(0.5)), Message(Gaussian(m=2.0, V=1.0))],
                                    Gaussian(m=2.0, V=1.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    1,
                                    [nothing, Message(MatrixDelta(0.5*diageye(2))), Message(MvGaussian(m=2.0*ones(2), V=1.0*diageye(2)))],
                                    MvGaussian(m=2.0*ones(2), V=3.0*diageye(2)))
        end

        context("Sum-product message towards variance/precision") do
            validateOutboundMessage(GaussianNode(),
                                    2,
                                    [Message(Delta(2.0)), nothing, Message(Delta(1.0))],
                                    InverseGamma(a=-0.5, b=0.5))
            validateOutboundMessage(GaussianNode(form=:precision),
                                    2,
                                    [Message(Delta(2.0)), nothing, Message(Delta(1.0))],
                                    Gamma(a=1.5, b=0.5))
        end
    end

    context("Sum-product: fixed mean") do
        context("Sum-product message towards out") do
            validateOutboundMessage(GaussianNode(m=2.0),
                                    2,
                                    [Message(Delta(0.5)), nothing],
                                    Gaussian(m=2.0, V=0.5))
            validateOutboundMessage(GaussianNode(m=2.0; form=:precision),
                                    2,
                                    [Message(Delta(0.5)), nothing],
                                    Gaussian(m=2.0, W=0.5))
        end

        context("Sum-product message towards variance/precision") do
            validateOutboundMessage(GaussianNode(m=2.0),
                                    1,
                                    [nothing, Message(Delta(1.0))],
                                    InverseGamma(a=-0.5, b=0.5))
            validateOutboundMessage(GaussianNode(m=2.0; form=:precision),
                                    1,
                                    [nothing, Message(Delta(1.0))],
                                    Gamma(a=1.5, b=0.5))
        end
    end

    context("Sum-product: fixed variance/precision") do
        context("Sum-product message towards out") do
            validateOutboundMessage(GaussianNode(V=0.5),
                                    2,
                                    [Message(Delta(2.0)), nothing],
                                    Gaussian(m=2.0, V=0.5))
        end

        context("Sum-product message towards mean") do
            validateOutboundMessage(GaussianNode(V=0.5),
                                    1,
                                    [nothing, Message(Delta(1.0))],
                                    Gaussian(m=1.0, V=0.5))
        end
    end

    context("Sum-product: all parameters fixed") do
            validateOutboundMessage(GaussianNode(m=2.0, V=0.5),
                                    1,
                                    [nothing],
                                    Gaussian(m=2.0, V=0.5))
    end

    context("Variational estimation") do
        context("GaussianNode with fixed variance should propagate a forward message to y") do
            validateOutboundMessage(GaussianNode(V=2.0),
                                    2,
                                    [Gaussian(m=3.0, V=1.0), nothing],
                                    Gaussian(m=3.0, V=2.0),
                                    ForneyLab.variationalRule!)
        end

        context("Naive variational implementation (mean field)") do
            context("GaussianNode should propagate a backward message to the mean") do
                # Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, Gamma(a=3.0, b=1.0), Gaussian(m=2.0, V=0.1)],
                                        Gaussian(m=2.0, W=3.0),
                                        ForneyLab.variationalRule!)
                # Mv Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, Wishart(V=[1.0 0.0; 0.0 2.0], nu=2.0), MvGaussian(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0])],
                                        MvGaussian(m=[1.0, 2.0], W=[2.0 0.0; 0.0 4.0]),
                                        ForneyLab.variationalRule!)
                # Variance
                validateOutboundMessage(GaussianNode(),
                                        1,
                                        [nothing, InverseGamma(a=3.0, b=1.0), Gaussian(m=2.0, V=0.1)],
                                        Gaussian(m=2.0, V=4.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a backward message to the variance, precision or log-precision") do
                # Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [Gaussian(m=4.0, W=2.0), nothing, Gaussian(m=2.0, V=0.1)],
                                        Gamma(a=1.5, b=2.3),
                                        ForneyLab.variationalRule!)
                # MvPrecision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [MvGaussian(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0]), nothing, MvGaussian(m=[3.0, 4.0], V=[1.0 0.0; 0.0 2.0])],
                                        Wishart(V=[0.25 -0.125; -0.125 0.1875], nu=4.0),
                                        ForneyLab.variationalRule!)
                # Variance
                validateOutboundMessage(GaussianNode(form=:variance),
                                        2,
                                        [Gaussian(m=4.0, V=1.0), nothing, Gaussian(m=2.0, V=0.1)],
                                        InverseGamma(a=-0.5, b=2.55),
                                        ForneyLab.variationalRule!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        1,
                                        [nothing, Gaussian(m=2.0, V=0.5)],
                                        Gamma(a=1.5, b=0.75),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=4.0; form=:variance),
                                        1,
                                        [nothing, Gaussian(m=2.0, V=0.1)],
                                        InverseGamma(a=-0.5, b=2.05),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a forward message") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [Gaussian(), Gamma(), nothing],
                                        Gaussian(m=0.0, W=1.0),
                                        ForneyLab.variationalRule!)
                # Mv Precision
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [MvGaussian(m=[1.0, 2.0], V=[1.0 0.0; 0.0 2.0]), Wishart(V=[1.0 0.0; 0.0 2.0], nu=2.0), nothing],
                                        MvGaussian(m=[1.0, 2.0], W=[2.0 0.0; 0.0 4.0]),
                                        ForneyLab.variationalRule!)
                # Variance
                validateOutboundMessage(GaussianNode(form=:variance),
                                        3,
                                        [Gaussian(), InverseGamma(a=2.0, b=1.0), nothing],
                                        Gaussian(m=0.0, V=1.0),
                                        ForneyLab.variationalRule!)
                # With fixed mean
                validateOutboundMessage(GaussianNode(m=1.0; form=:precision),
                                        2,
                                        [Gamma(a=1.5, b=0.5), nothing],
                                        Gaussian(m=1.0, W=3.0),
                                        ForneyLab.variationalRule!)
                validateOutboundMessage(GaussianNode(m=1.0; form=:variance),
                                        2,
                                        [InverseGamma(a=2.0, b=1.0), nothing],
                                        Gaussian(m=1.0, V=1.0),
                                        ForneyLab.variationalRule!)
            end
        end

        context("Structured variational implementation") do
            context("GaussianNode should propagate a backward message to the mean") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        1,
                                        [nothing, Message(Gamma(a=3.0, b=1.0)), Gaussian(m=2.0, W=10.0)],
                                        StudentsT(m=2.0, lambda=60.0/21.0, nu=6.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a backward message to the precision") do
                validateOutboundMessage(GaussianNode(form=:precision),
                                        2,
                                        [Message(Gaussian(m=4.0, W=2.0)), nothing, Gaussian(m=2.0, W=10.0)],
                                        Gamma(a=1.5, b=41.0/20.0),
                                        ForneyLab.variationalRule!)
            end

            context("GaussianNode should propagate a forward message") do
                node_marg = NormalGamma(m=4.0, beta=2.0, a=3.0, b=1.0)
                validateOutboundMessage(GaussianNode(form=:precision),
                                        3,
                                        [node_marg, node_marg, nothing],
                                        Gaussian(m=4.0, W=3.0),
                                        ForneyLab.variationalRule!)
            end
        end
    end
end
