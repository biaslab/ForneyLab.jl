#####################
# Unit tests
#####################

facts("GaussianMixtureNode unit tests") do
    context("GaussianmixtureNode should initialize a GaussianMixtureNode with 4 interfaces") do
        FactorGraph()
        GaussianMixtureNode(id=:node)
        @fact length(n(:node).interfaces) --> 4
        @fact n(:node).i[:m] --> n(:node).interfaces[1]
        @fact n(:node).i[:w] --> n(:node).interfaces[2]
        @fact n(:node).i[:x] --> n(:node).interfaces[3]
        @fact n(:node).i[:z] --> n(:node).interfaces[4]
    end

    context("GaussianMixtureNode should pass variational messages") do
      # message to m
      validateOutboundMessage(GaussianMixtureNode(),
                              1,
                              [nothing,Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.)]), Gaussian(m=8.5,V=0.5), Bernoulli(0.2)],
                              Partitioned([Gaussian(m=8.5,V=10.0),Gaussian(m=8.5,V=0.625)]),
                              ForneyLab.variationalRule!)

      # message to lambda
      validateOutboundMessage(GaussianMixtureNode(),
                              2,
                              [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.)]),nothing, Gaussian(m=8.5,V=0.5), Bernoulli(0.2)],
                              Partitioned([Gamma(a=1.1,b=1.475), Gamma(a=1.4,b=2.3000000000000003)]),
                              ForneyLab.variationalRule!)

      # message to x
      validateOutboundMessage(GaussianMixtureNode(),
                              3,
                              [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.)]),Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.)]), nothing, Bernoulli(0.2)],
                              Gaussian(m=9.705882352941178,V=0.5882352941176471),
                              ForneyLab.variationalRule!)

        # message to z
        validateOutboundMessage(GaussianMixtureNode(),
                                4,
                                [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.)]),Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.)]), Gaussian(m=8.5,V=0.5), nothing],
                                Bernoulli(0.7713458788198755),
                                ForneyLab.variationalRule!)

        #message to m univariate multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                1,
                                [nothing,Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.), Gamma(a=3.,b=1.)]), Gaussian(m=8.5,V=0.5), Categorical([0.1,0.1,0.8])],
                                Partitioned([Gaussian(m=8.5,V=20.),Gaussian(m=8.5,V=5.),Gaussian(m=8.5,V=0.41666666666666663)]),
                                ForneyLab.variationalRule!)

        #message to lambda univariate multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                2,
                                [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.),Gaussian(m=12.,V=3.)]),nothing, Gaussian(m=8.5,V=0.5), Categorical([0.1,0.1,0.8])],
                                Partitioned([Gamma(a=1.05,b=0.7375), Gamma(a=1.05,b=0.28750000000000003), Gamma(a=1.4,b=6.300000000000001)]),
                                ForneyLab.variationalRule!)

        #message to x univariate multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                3,
                                [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.),Gaussian(m=12.,V=3.)]),Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.), Gamma(a=3.,b=1.)]), nothing, Categorical([0.1,0.1,0.8])],
                                Gaussian(m=11.716981132075471,V=0.3773584905660377),
                                ForneyLab.variationalRule!)

        #message to z univariate multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                4,
                                [Partitioned([Gaussian(m=5.,V=2.),Gaussian(m=10.,V=3.),Gaussian(m=12.,V=3.)]),Partitioned([Gamma(a=1.,b=2.), Gamma(a=2.,b=1.), Gamma(a=3.,b=1.)]), Gaussian(m=8.5,V=0.5), nothing],
                                Categorical([0.7713458749115749,0.2286541200215665,5.0668586342559045e-9]),
                                ForneyLab.variationalRule!)


        #message to m multivariate with two clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                1,
                                [nothing,Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.)]), MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), Bernoulli(0.2)],
                                Partitioned([MvGaussian(m=[8.5,1.],V=[1.6666666666666665 0.0; 0.0 1.6666666666666665]),MvGaussian(m=[8.5,1.],V=[0.20833333333333331 0.0; 0.0 0.20833333333333331])]),
                                ForneyLab.variationalRule!)

        #message to lambda multivariate with two clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                2,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.])]),nothing, MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), Bernoulli(0.2)],
                                Partitioned([Wishart(V=[0.3837953091684434 -0.2985074626865672;-0.2985074626865672 1.3432835820895523],nu=3.2), Wishart(V=[0.28225806451612895 -0.12096774193548387;-0.12096774193548387 0.7661290322580646],nu=3.8)]),
                                ForneyLab.variationalRule!)

        #message to x multivariate with two clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                3,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.])]),Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.)]), nothing, Bernoulli(0.2)],
                                MvGaussian(xi=[51.0,7.2],V=[0.18518518518518517 0;0 0.18518518518518517]),
                                ForneyLab.variationalRule!)

        #message to z multivariate with two clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                4,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.])]),Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.)]), MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), nothing],
                                Bernoulli(9.462663664175011e-6),
                                ForneyLab.variationalRule!)

        #message to m multivariate with multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                1,
                                [nothing,Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.),Wishart(V=[2. 0.;0. 2.],nu=4.)]), MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), Categorical([0.1, 0.1, 0.8])],
                                Partitioned([MvGaussian(m=[8.5,1.],V=[3.333333333333333 0.0; 0.0 3.333333333333333]),MvGaussian(m=[8.5,1.],V=[1.6666666666666665 0.0; 0.0 1.6666666666666665]),MvGaussian(m=[8.5,1.],V=[0.15625 0.0; 0.0 0.15625])]),
                                ForneyLab.variationalRule!)

        #message to lambda multivariate with multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                2,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.]), MvGaussian(m=[12.,1.],V=[2. 0.0; 0.0 1.])]),nothing, MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), Categorical([0.1, 0.1, 0.8])],
                                Partitioned([Wishart(V=[0.7675906183368868 -0.5970149253731344;-0.5970149253731339 2.6865671641791047],nu=3.1), Wishart(V=[2.2580645161290316 -0.967741935483871;-0.9677419354838712 6.129032258064517],nu=3.1),Wishart(V=[0.0847457627118644 0.;0. 0.8333333333333333],nu=3.8)]),
                                ForneyLab.variationalRule!)

        #message to lambda multivariate with multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                3,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.]), MvGaussian(m=[12.,1.],V=[2. 0.0; 0.0 1.])]),Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.),Wishart(V=[2. 0.;0. 2.],nu=4.)]), nothing, Categorical([0.1, 0.1, 0.8])],
                                MvGaussian(xi = [84.30000000000001,7.300000000000001],V=[0.136986301369863 0;0 0.136986301369863]),
                                ForneyLab.variationalRule!)

        #message to z multivariate with multiple clusters
        validateOutboundMessage(GaussianMixtureNode(),
                                4,
                                [Partitioned([MvGaussian(m=[5.,0.],V=[3 0.0; 0.0 3.]),MvGaussian(m=[10.,1.5],V=[2. 0.0; 0.0 1.]), MvGaussian(m=[12.,1.],V=[2. 0.0; 0.0 1.])]),Partitioned([Wishart(V=[1. 0.;0. 1.],nu=3.), Wishart(V=[2. 0.;0. 2.],nu=3.),Wishart(V=[2. 0.;0. 2.],nu=4.)]),  MvGaussian(m=[8.5,1.],V=[0.5 0;0 0.5]), nothing],
                                Categorical([9.462663664175011e-6, 0.9999905373363358, 2.8624914934988436e-20]),
                                ForneyLab.variationalRule!)

    end
end

# #####################
# # Integration tests
# #####################
# facts("CategoricalNode integration tests") do
#     context("Variational message passing should estimate Dirichlet distribution parameters") do
#         data = [Categorical([0.1,0.1,0.8]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1])]
#         K = length(data)
#
#         FactorGraph()
#
#         for k = 1:K
#             EqualityNode(id=:eq_*k)
#             CategoricalNode(id=:bern_*k)
#             TerminalNode(data[k], id=:obs_*k)
#
#             Edge(n(:eq_*k).i[3], n(:bern_*k).i[:pi], id=:x_*k)
#             Edge(n(:bern_*k).i[:z], n(:obs_*k).i[:out], id=:c_*k)
#
#             if k > 1
#                 Edge(n(:eq_*(k-1)).i[2], n(:eq_*k).i[1])
#             end
#         end
#
#         TerminalNode(Dirichlet([1.,1.,1.]), id=:prior)
#         TerminalNode(vague(Dirichlet{3}), id=:term)
#
#         Edge(n(:prior).i[:out], n(:eq_1).i[1])
#         Edge(n(:eq_*K).i[2], n(:term).i[:out])
#
#         buff = attachWriteBuffer(n(:term).i[:out].partner)
#
#         RecognitionFactorization()
#
#         factor(eg(:x_1))
#         for k = 1:K
#             factor(eg(:c_*k))
#
#             initialize(eg(:c_*k), vague(Categorical{3}))
#             initialize(eg(:x_*k), vague(Dirichlet{3}))
#         end
#
#         algo = VariationalBayes(n_iterations=50)
#         run(algo)
#
#         @fact buff[end].alpha --> [6.,1.,1.]
#     end
# end
