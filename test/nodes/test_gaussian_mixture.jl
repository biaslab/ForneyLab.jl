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
                                Bernoulli(0.0006053637068648643),
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
                                Categorical([0.0006049974633578021, 0.9987900050732843, 0.0006049974633578021]),
                                ForneyLab.variationalRule!)

    end
end

# #####################
# # Integration tests
# #####################
facts("Gaussian Mixture integration tests") do
    context("Variational message passing should estimate messages, univariate two clusters") do
      data = [4.71772,20.6781,22.5677,19.0028,1.32842,2.80178,3.13219,3.64153] # data sampled from two normal distributions with mean 20.0 and 3.0
      N    = length(data)
      n_its = 50

      FactorGraph()

      for k=1:(N)
      GaussianMixtureNode(id=:gm*k)
      BernoulliNode(id=:ber*k)
      Edge(n(:ber*k).i[:out],n(:gm*k).i[:z],id=:z_e*k)
      TerminalNode(Gaussian(m=data[k],V=0.5),id=:y*k)
      Edge(n(:y*k).i[:out],n(:gm*k).i[:x],id=:y_e*k)
      end

      PriorNode(Partitioned([Gaussian(m=23.,V=10.),Gaussian(m=9.0,V=10.)]),id=:m1_start)
      PriorNode(Partitioned([ForneyLab.Gamma(a=10.,b=1.),ForneyLab.Gamma(a=10.,b=1.)]),id=:w1_start)
      PriorNode(ForneyLab.Beta(a=2.,b=2.),id=:pi_start)
      for k = 1:(N-1)
      EqualityNode(id=:m1_eq*k)
      EqualityNode(id=:w1_eq*k)
      EqualityNode(id=:pi_eq*k)
      end

      Edge(n(:gm*1).i[:m],n(:m1_eq*1).i[:3],id=:m_e*1)
      Edge(n(:gm*1).i[:w],n(:w1_eq*1).i[:3],id=:w_e*1)
      Edge(n(:ber*1).i[:in],n(:pi_eq*1).i[:3],id=:pi_e*1)

      for k = 2:N-1
        Edge(n(:gm*k).i[:m],n(:m1_eq*k).i[:3],id=:m_e*k)
        Edge(n(:m1_eq*(k-1)).i[:2],n(:m1_eq*k).i[:1])

        Edge(n(:gm*k).i[:w],n(:w1_eq*k).i[:3],id=:w_e*k)
        Edge(n(:w1_eq*(k-1)).i[:2],n(:w1_eq*k).i[:1])

        Edge(n(:ber*k).i[:in],n(:pi_eq*k).i[:3],id=:pi_e*k)
        Edge(n(:pi_eq*(k-1)).i[:2],n(:pi_eq*k).i[:1])
      end
      Edge(n(:gm*(N)).i[:m],n(:m1_eq*(N-1)).i[:2],id=:m_e*(N))
      Edge(n(:m1_start).i[:out],n(:m1_eq*1).i[:1],id=:m_e_s*1)

      Edge(n(:gm*(N)).i[:w],n(:w1_eq*(N-1)).i[:2],id=:w_e*(N))
      Edge(n(:w1_start).i[:out],n(:w1_eq*1).i[:1],id=:w_e_s*1)


      Edge(n(:ber*(N)).i[:in],n(:pi_eq*(N-1)).i[:2],id=:pi_e*(N))
      Edge(n(:pi_start).i[:out],n(:pi_eq*1).i[:1],id=:pi_e_s*1)

      #Attach write buffers
      m1_est = attachWriteBuffer(n(:m1_start).i[:out].partner)
      w1_est = attachWriteBuffer(n(:w1_start).i[:out].partner)
      pi_est = attachWriteBuffer(n(:pi_start).i[:out].partner)
      z_est = attachWriteBuffer(n(:ber*1).i[:out].partner)

      # Specify the variational algorithm for n_its vmp iterations
      rf = RecognitionFactorization()
      factorizeMeanField(rf)
      for k=1:(N)
      initialize(eg(:pi_e*k), ForneyLab.Beta(a=2.,b=2.))
      initialize(eg(:m_e*k), Partitioned([Gaussian(m=23.,V=10.),Gaussian(m=9.0,V=10.)]))
      initialize(eg(:w_e*k), Partitioned([ForneyLab.Gamma(a=10.,b=1.),ForneyLab.Gamma(a=10.,b=1.)]))
      initialize(eg(:y_e*k),Gaussian(m=data[k],V=0.5))
      initialize(eg(:z_e*k),ForneyLab.Bernoulli(0.5))
      end

      # Specify the variational algorithm for n_its vmp iterations
      msg_types = Dict{Interface,DataType}(n(:gm*i).i[:z] => ForneyLab.Bernoulli for i=1:(N))

      algo = VariationalBayes(n_iterations=n_its, message_types=msg_types)

      run(algo)
      ensureParameters!(m1_est[end].factors[1], (:m, :V))
      ensureParameters!(m1_est[end].factors[2], (:m, :V))

      @fact pi_est[end].a/(pi_est[end].a+pi_est[end].b) --> less_than(0.4)
      @fact pi_est[end].a/(pi_est[end].a+pi_est[end].b) --> greater_than(0.399)
      @fact m1_est[end].factors[1].m --> greater_than(20.786)
      @fact m1_est[end].factors[1].m --> less_than(20.79)
      @fact m1_est[end].factors[2].m --> greater_than(3.18)
      @fact m1_est[end].factors[2].m --> less_than(3.19)
      @fact m1_est[end].factors[1].V --> roughly(0.03779450341624337)
      @fact m1_est[end].factors[2].V --> roughly(0.022123670531282245)
      @fact z_est[1].p--> less_than(0.000001)
      @fact w1_est[end].factors[1].a/w1_est[end].factors[1].b --> roughly(8.226109610890315)
      @fact w1_est[end].factors[2].a/w1_est[end].factors[2].b --> roughly(9.144842051081852)
     end

     context("Variational message passing should estimate messages, univariate multiple clusters") do
       data =[-10.6962; 3.94545; 3.62594; -10.3669; 26.0386; 2.91566; -8.34216; 26.9278; 3.49497; 3.24912; 21.2678; 17.1254; -11.48; 19.5314; -9.80081] # data sampled from two normal distributions with mean 20.0 and 3.0
       N    = length(data)
       n_its = 50

       FactorGraph()

       for k=1:(N)
       GaussianMixtureNode(id=:gm*k)
       CategoricalNode(id=:ber*k)
       Edge(n(:ber*k).i[:z],n(:gm*k).i[:z],id=:z_e*k)
       TerminalNode(Gaussian(m=data[k],V=0.5),id=:y*k)
       Edge(n(:y*k).i[:out],n(:gm*k).i[:x],id=:y_e*k)
       end

       PriorNode(Partitioned([Gaussian(m=23.,V=10.),Gaussian(m=9.,V=4.),Gaussian(m=-2.,V=3.)]),id=:m1_start)
       PriorNode(Partitioned([ForneyLab.Gamma(a=1., b=1.),ForneyLab.Gamma(a=1., b=1.),ForneyLab.Gamma(a=1., b=1.)]),id=:w1_start)
       PriorNode(ForneyLab.Dirichlet([2.,2.,2.]),id=:pi_start)
       for k = 1:(N-1)
       EqualityNode(id=:m1_eq*k)
       EqualityNode(id=:w1_eq*k)
       EqualityNode(id=:pi_eq*k)
       end

       Edge(n(:gm*1).i[:m],n(:m1_eq*1).i[:3],id=:m_e*1)
       Edge(n(:gm*1).i[:w],n(:w1_eq*1).i[:3],id=:w_e*1)
       Edge(n(:ber*1).i[:pi],n(:pi_eq*1).i[:3],id=:pi_e*1)

       for k = 2:N-1
         Edge(n(:gm*k).i[:m],n(:m1_eq*k).i[:3],id=:m_e*k)
         Edge(n(:m1_eq*(k-1)).i[:2],n(:m1_eq*k).i[:1])

         Edge(n(:gm*k).i[:w],n(:w1_eq*k).i[:3],id=:w_e*k)
         Edge(n(:w1_eq*(k-1)).i[:2],n(:w1_eq*k).i[:1])

         Edge(n(:ber*k).i[:pi],n(:pi_eq*k).i[:3],id=:pi_e*k)
         Edge(n(:pi_eq*(k-1)).i[:2],n(:pi_eq*k).i[:1])
       end
       Edge(n(:gm*(N)).i[:m],n(:m1_eq*(N-1)).i[:2],id=:m_e*(N))
       Edge(n(:m1_start).i[:out],n(:m1_eq*1).i[:1],id=:m_e_s*1)

       Edge(n(:gm*(N)).i[:w],n(:w1_eq*(N-1)).i[:2],id=:w_e*(N))
       Edge(n(:w1_start).i[:out],n(:w1_eq*1).i[:1],id=:w_e_s*1)


       Edge(n(:ber*(N)).i[:pi],n(:pi_eq*(N-1)).i[:2],id=:pi_e*(N))
       Edge(n(:pi_start).i[:out],n(:pi_eq*1).i[:1],id=:pi_e_s*1)

       #Attach write buffers
       m1_est = attachWriteBuffer(n(:m1_start).i[:out].partner)
       w1_est = attachWriteBuffer(n(:w1_start).i[:out].partner)
       pi_est = attachWriteBuffer(n(:pi_start).i[:out].partner)
       z_est = attachWriteBuffer(n(:ber*1).i[:z].partner)

       # Specify the variational algorithm for n_its vmp iterations
       rf = RecognitionFactorization()
       factorizeMeanField(rf)
       for k=1:(N)
       initialize(eg(:pi_e*k), ForneyLab.Dirichlet([2.,2.,2.]))
       initialize(eg(:m_e*k), Partitioned([Gaussian(m=23.,V=10.),Gaussian(m=9.,V=4.),Gaussian(m=-2.,V=3.)]))
       initialize(eg(:w_e*k), Partitioned([ForneyLab.Gamma(a=1., b=1.),ForneyLab.Gamma(a=1., b=1.),ForneyLab.Gamma(a=1., b=1.)]))
       initialize(eg(:y_e*k),Gaussian(m=data[k],V=0.5))
       initialize(eg(:z_e*k),ForneyLab.Categorical([0.25,0.3,0.45]))
       end

       # Specify the variational algorithm for n_its vmp iterations
       msg_types = Dict{Interface,DataType}(n(:gm*i).i[:z] => ForneyLab.Categorical{3} for i=1:(N))

       algo = VariationalBayes(n_iterations=n_its, message_types=msg_types)

       run(algo)
       ensureParameters!(m1_est[end].factors[1], (:m, :V))
       ensureParameters!(m1_est[end].factors[2], (:m, :V))
       ensureParameters!(m1_est[end].factors[3], (:m, :V))

       @fact pi_est[end].alpha[1]/sum(pi_est[end].alpha) --> greater_than(0.3333)
       @fact pi_est[end].alpha[1]/sum(pi_est[end].alpha) --> less_than(0.334)
       @fact m1_est[end].factors[1].m -->greater_than(22.18)
       @fact m1_est[end].factors[1].m -->less_than(22.19)
       @fact m1_est[end].factors[2].m --> greater_than(3.57)
       @fact m1_est[end].factors[2].m --> less_than(3.58)
       @fact m1_est[end].factors[3].m --> greater_than(-9.55)
       @fact m1_est[end].factors[3].m --> less_than(-9.44)
       @fact m1_est[end].factors[1].V --> greater_than(2.24)
       @fact m1_est[end].factors[1].V --> less_than(2.25)
       @fact m1_est[end].factors[2].V --> greater_than(0.13)
       @fact m1_est[end].factors[2].V --> less_than(0.14)
       @fact m1_est[end].factors[3].V --> greater_than(4.00)
       @fact m1_est[end].factors[3].V --> less_than(4.05)
       @fact z_est[1].p[1]--> less_than(0.00001)
       @fact w1_est[end].factors[1].a/w1_est[end].factors[1].b --> greater_than(0.09)
       @fact w1_est[end].factors[1].a/w1_est[end].factors[1].b --> less_than(0.095)
       @fact w1_est[end].factors[2].a/w1_est[end].factors[2].b --> greater_than(2.93)
       @fact w1_est[end].factors[2].a/w1_est[end].factors[2].b --> less_than(2.94)
       @fact w1_est[end].factors[3].a/w1_est[end].factors[3].b --> greater_than(0.048)
       @fact w1_est[end].factors[3].a/w1_est[end].factors[3].b --> less_than(0.049)
      end

      context("Variational message passing should estimate messages, multivariate 2 clusters") do
        data = [20.2366 3.57175 4.80739; 2.27229 7.57345 20.2518; 3.85124 6.3984 20.2731; 1.73924 10.818 20.1451; 20.1333 6.29229 6.68207] # data sampled from two normal distributions with mean 20.0 and 3.0
        N    = size(data)[1]
        n_its = 50

        FactorGraph()

        for k=1:(N)
        GaussianMixtureNode(id=:gm*k)
        BernoulliNode(id=:ber*k)
        Edge(n(:ber*k).i[:out],n(:gm*k).i[:z],id=:z_e*k)
        TerminalNode(MvGaussian(m=reshape(data[k,:],size(data,2)),V=[0.5 0. 0.; 0. 0.5 0.;0. 0. 0.5]),id=:y*k)
        Edge(n(:y*k).i[:out],n(:gm*k).i[:x],id=:y_e*k)
        end

        PriorNode(Partitioned([MvGaussian(m=[23.,10.0, 1.0],V=[10. 1. 1.;1. 10. 1.; 1. 1. 10.0]),MvGaussian(m=[9.0, 10.0, 15.0],V=[10. 1. 1.; 1. 10. 1.; 1. 1. 10.0])]),id=:m1_start)
        PriorNode(Partitioned([ForneyLab.Wishart(nu=3., V=eye(3)+0.05*ones(3,3)),ForneyLab.Wishart(nu=3., V=eye(3)+0.05*ones(3,3))]),id=:w1_start)
        PriorNode(ForneyLab.Beta(a=2.,b=2.),id=:pi_start)
        for k = 1:(N-1)
        EqualityNode(id=:m1_eq*k)
        EqualityNode(id=:w1_eq*k)
        EqualityNode(id=:pi_eq*k)
        end

        Edge(n(:gm*1).i[:m],n(:m1_eq*1).i[:3],id=:m_e*1)
        Edge(n(:gm*1).i[:w],n(:w1_eq*1).i[:3],id=:w_e*1)
        Edge(n(:ber*1).i[:in],n(:pi_eq*1).i[:3],id=:pi_e*1)

        for k = 2:N-1
          Edge(n(:gm*k).i[:m],n(:m1_eq*k).i[:3],id=:m_e*k)
          Edge(n(:m1_eq*(k-1)).i[:2],n(:m1_eq*k).i[:1])

          Edge(n(:gm*k).i[:w],n(:w1_eq*k).i[:3],id=:w_e*k)
          Edge(n(:w1_eq*(k-1)).i[:2],n(:w1_eq*k).i[:1])

          Edge(n(:ber*k).i[:in],n(:pi_eq*k).i[:3],id=:pi_e*k)
          Edge(n(:pi_eq*(k-1)).i[:2],n(:pi_eq*k).i[:1])
        end
        Edge(n(:gm*(N)).i[:m],n(:m1_eq*(N-1)).i[:2],id=:m_e*(N))
        Edge(n(:m1_start).i[:out],n(:m1_eq*1).i[:1],id=:m_e_s*1)

        Edge(n(:gm*(N)).i[:w],n(:w1_eq*(N-1)).i[:2],id=:w_e*(N))
        Edge(n(:w1_start).i[:out],n(:w1_eq*1).i[:1],id=:w_e_s*1)


        Edge(n(:ber*(N)).i[:in],n(:pi_eq*(N-1)).i[:2],id=:pi_e*(N))
        Edge(n(:pi_start).i[:out],n(:pi_eq*1).i[:1],id=:pi_e_s*1)

        #Attach write buffers
        m1_est = attachWriteBuffer(n(:m1_start).i[:out].partner)
        w1_est = attachWriteBuffer(n(:w1_start).i[:out].partner)
        pi_est = attachWriteBuffer(n(:pi_start).i[:out].partner)
        z_est = attachWriteBuffer(n(:ber*1).i[:out].partner)

        # Specify the variational algorithm for n_its vmp iterations
        rf = RecognitionFactorization()
        factorizeMeanField(rf)
        for k=1:(N)
          initialize(eg(:pi_e*k), ForneyLab.Beta(a=2.,b=2.))
          initialize(eg(:m_e*k), Partitioned([MvGaussian(m=[23.,10.0, 1.0],V=[10. 0.1 0.1;0.1 10. 0.1; 0.1 0.1 10.0]),MvGaussian(m=[9.0, 10.0, 15.0],V=[10. 0.1 0.1; 0.1 10. 0.1; 0.1 0.1 10.0])]))
          initialize(eg(:w_e*k), Partitioned([ForneyLab.Wishart(nu=3.,  V=eye(3)/3.),ForneyLab.Wishart(nu=3.,  V=eye(3)/3.)]))
          initialize(eg(:y_e*k),MvGaussian(m=reshape(data[k,:],size(data,2)),V=[0.5 0.1 0.1; 0.1 0.5 0.1;0.1 0.1 0.5]))
          initialize(eg(:z_e*k),ForneyLab.Bernoulli(0.5))
        end

        # Specify the variational algorithm for n_its vmp iterations
        msg_types = Dict{Interface,DataType}(n(:gm*i).i[:z] => ForneyLab.Bernoulli for i=1:(N))

        algo = VariationalBayes(n_iterations=n_its, message_types=msg_types)

        run(algo)
        ensureParameters!(m1_est[end].factors[1], (:m, :V))
        ensureParameters!(m1_est[end].factors[2], (:m, :V))


        @fact pi_est[end].a/sum(pi_est[end].a+pi_est[end].b) --> greater_than(0.42)
        @fact pi_est[end].a/sum(pi_est[end].a+pi_est[end].b) --> less_than(0.43)
        @fact m1_est[end].factors[1].m[1] --> greater_than(20.24)
        @fact m1_est[end].factors[1].m[1] -->less_than(20.26)
        @fact m1_est[end].factors[2].m[2] --> greater_than(8.27)
        @fact m1_est[end].factors[2].m[2] --> less_than(8.30)
        @fact m1_est[end].factors[1].V[1] --> greater_than(0.17)
        @fact m1_est[end].factors[1].V[1] --> less_than(0.18)
        @fact m1_est[end].factors[2].V[3] --> greater_than(0.0006)
        @fact m1_est[end].factors[2].V[3] --> less_than(0.00065)
        @fact z_est[1].p[1]--> greater_than(0.999)
        @fact w1_est[end].factors[1].nu*w1_est[end].factors[1].V[1] --> greater_than(8.05)
        @fact w1_est[end].factors[1].nu*w1_est[end].factors[1].V[1] --> less_than(8.1)
        @fact w1_est[end].factors[2].nu*w1_est[end].factors[2].V[2] --> greater_than(2.32)
        @fact w1_est[end].factors[2].nu*w1_est[end].factors[2].V[2] --> less_than(2.36)
       end

       context("Variational message passing should estimate messages, multivariate multiple clusters") do
         data =[-3.31872 10.0357 30.1384; 3.47228 11.9783 20.379; 20.1562 6.37659 7.37432; -5.0689 9.35502 29.8067; 0.743716 12.8442 18.5372; 19.1677 5.56954 6.38645; 20.4842 4.92519 7.08919; 19.9359 5.90297 7.87806; 2.44867 12.2538 20.8173; 4.78051 10.6962 20.624; -3.38729 10.8782 29.9854; 2.69968 11.2883 20.2603; 19.8764 5.0597 4.84675; 4.86845 10.6923 20.4384; -5.0638 9.47881 29.3086; -3.7251 10.193 31.9495; 21.1899 5.36709 4.59334; 3.44318 9.08621 18.4571; 20.8348 4.66949 5.76017; -4.75307 12.4824 29.4761; 19.5512 5.81893 6.83283; -5.13072 10.9483 30.256; -5.32767 9.66011 30.9059; 2.74185 6.68038 18.1188; -5.17731 9.17734 30.7458; -4.68014 9.00945 29.6801; 23.8108 5.44229 5.95253; 23.0177 4.28682 5.13516; 3.07371 9.22484 20.8193; 2.51827 7.35499 18.3339] # data sampled from two normal distributions with mean 20.0 and 3.0
         N    = size(data)[1]
         n_its = 10

         FactorGraph()

         for k=1:(N)
         GaussianMixtureNode(id=:gm*k)
         CategoricalNode(id=:ber*k)
         Edge(n(:ber*k).i[:z],n(:gm*k).i[:z],id=:z_e*k)
         TerminalNode(MvGaussian(m=reshape(data[k,:],size(data,2)),V=0.5*eye(3)),id=:y*k)
         Edge(n(:y*k).i[:out],n(:gm*k).i[:x],id=:y_e*k)
         end

         PriorNode(Partitioned([MvGaussian(m=[23.,7.0, 5.0],V=[10. 1. 1.;1. 10. 1.; 1. 1. 10.0]),MvGaussian(m=[5.0, 10.0, 18.0],V=[10. 1. 1.; 1. 10. 1.; 1. 1. 10.0]),MvGaussian(m=[-3.0, 10.0, 28.0],V=[10. 1. 1.; 1. 10. 1.; 1. 1. 10.0])]),id=:m1_start)
         PriorNode(Partitioned([ForneyLab.Wishart(nu=3., V=eye(3)),ForneyLab.Wishart(nu=3., V=eye(3)),ForneyLab.Wishart(nu=3., V=eye(3))]),id=:w1_start)
         PriorNode(ForneyLab.Dirichlet([2.,2.,2.]),id=:pi_start)
         for k = 1:(N-1)
         EqualityNode(id=:m1_eq*k)
         EqualityNode(id=:w1_eq*k)
         EqualityNode(id=:pi_eq*k)
         end

         Edge(n(:gm*1).i[:m],n(:m1_eq*1).i[:3],id=:m_e*1)
         Edge(n(:gm*1).i[:w],n(:w1_eq*1).i[:3],id=:w_e*1)
         Edge(n(:ber*1).i[:pi],n(:pi_eq*1).i[:3],id=:pi_e*1)

         for k = 2:N-1
           Edge(n(:gm*k).i[:m],n(:m1_eq*k).i[:3],id=:m_e*k)
           Edge(n(:m1_eq*(k-1)).i[:2],n(:m1_eq*k).i[:1])

           Edge(n(:gm*k).i[:w],n(:w1_eq*k).i[:3],id=:w_e*k)
           Edge(n(:w1_eq*(k-1)).i[:2],n(:w1_eq*k).i[:1])

           Edge(n(:ber*k).i[:pi],n(:pi_eq*k).i[:3],id=:pi_e*k)
           Edge(n(:pi_eq*(k-1)).i[:2],n(:pi_eq*k).i[:1])
         end
         Edge(n(:gm*(N)).i[:m],n(:m1_eq*(N-1)).i[:2],id=:m_e*(N))
         Edge(n(:m1_start).i[:out],n(:m1_eq*1).i[:1],id=:m_e_s*1)

         Edge(n(:gm*(N)).i[:w],n(:w1_eq*(N-1)).i[:2],id=:w_e*(N))
         Edge(n(:w1_start).i[:out],n(:w1_eq*1).i[:1],id=:w_e_s*1)


         Edge(n(:ber*(N)).i[:pi],n(:pi_eq*(N-1)).i[:2],id=:pi_e*(N))
         Edge(n(:pi_start).i[:out],n(:pi_eq*1).i[:1],id=:pi_e_s*1)

         #Attach write buffers
         m1_est = attachWriteBuffer(n(:m1_start).i[:out].partner)
         w1_est = attachWriteBuffer(n(:w1_start).i[:out].partner)
         pi_est = attachWriteBuffer(n(:pi_start).i[:out].partner)
         z_est = attachWriteBuffer(n(:ber*1).i[:z].partner)

         # Specify the variational algorithm for n_its vmp iterations
         rf = RecognitionFactorization()
         factorizeMeanField(rf)
         for k=1:(N)
           initialize(eg(:pi_e*k), ForneyLab.Dirichlet([2.,2.,2.]))
           initialize(eg(:m_e*k), Partitioned([MvGaussian(m=[23.,10.0, 1.0],V=[10. 0.1 0.1;0.1 10. 0.1; 0.1 0.1 10.0]),MvGaussian(m=[9.0, 10.0, 15.0],V=[10. 0.1 0.1; 0.1 10. 0.1; 0.1 0.1 10.0]),MvGaussian(m=[0.0, 10.0, 36.0],V=[10. 0.1 0.1; 0.1 10. 0.1; 0.1 0.1 10.0])]))
           initialize(eg(:w_e*k), Partitioned([ForneyLab.Wishart(nu=3., V=eye(3)),ForneyLab.Wishart(nu=3., V=eye(3)),ForneyLab.Wishart(nu=3., V=eye(3))]))
           initialize(eg(:y_e*k),MvGaussian(m=reshape(data[k,:],size(data,2)),V=[0.5 0.1 0.1; 0.1 0.5 0.1;0.1 0.1 0.5]))
           initialize(eg(:z_e*k),ForneyLab.Categorical([1/3, 1/3, 1/3]))
         end

         # Specify the variational algorithm for n_its vmp iterations
         msg_types = Dict{Interface,DataType}(n(:gm*i).i[:z] => ForneyLab.Categorical{3} for i=1:(N))

         algo = VariationalBayes(n_iterations=n_its, message_types=msg_types)

         run(algo)
         ensureParameters!(m1_est[end].factors[1], (:m, :V))
         ensureParameters!(m1_est[end].factors[2], (:m, :V))
         ensureParameters!(m1_est[end].factors[3], (:m, :V))

         @fact pi_est[end].alpha[1]/sum(pi_est[end].alpha) --> greater_than(0.3333)
         @fact pi_est[end].alpha[1]/sum(pi_est[end].alpha) --> less_than(0.334)
         @fact m1_est[end].factors[1].m[1] -->greater_than(20.81)
         @fact m1_est[end].factors[1].m[1] -->less_than(20.85)
         @fact m1_est[end].factors[2].m[2] --> greater_than(10.21)
         @fact m1_est[end].factors[2].m[2] --> less_than(10.22)
         @fact m1_est[end].factors[3].m[3] --> greater_than(30.1)
         @fact m1_est[end].factors[3].m[3] --> less_than(30.2)
         @fact m1_est[end].factors[1].V[1] --> greater_than(0.132)
         @fact m1_est[end].factors[1].V[1] --> less_than(0.133)
         @fact m1_est[end].factors[2].V[2] --> greater_than(-0.023)
         @fact m1_est[end].factors[2].V[2] --> less_than(-0.022)
         @fact m1_est[end].factors[3].V[3] --> greater_than(-0.101)
         @fact m1_est[end].factors[3].V[3] --> less_than(-0.100)
         @fact z_est[1].p[1]--> less_than(0.00001)
         @fact w1_est[end].factors[1].nu*w1_est[end].factors[1].V[1] --> greater_than(1.30)
         @fact w1_est[end].factors[1].nu*w1_est[end].factors[1].V[1] --> less_than(1.31)
         @fact w1_est[end].factors[2].nu*w1_est[end].factors[2].V[2] --> greater_than(0.54)
         @fact w1_est[end].factors[2].nu*w1_est[end].factors[2].V[2] --> less_than(0.55)
         @fact w1_est[end].factors[3].nu*w1_est[end].factors[3].V[3]  --> greater_than(0.95)
         @fact w1_est[end].factors[3].nu*w1_est[end].factors[3].V[3]  --> less_than(0.96)
        end

 end
