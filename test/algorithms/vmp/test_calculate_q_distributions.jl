    # context("Q distribution calculation for naively factorized GaussianNode") do
    #     (node, edges) = initializeGaussianNode(y_type=Float64)
    #     scheme = InferenceScheme()
    #     factorize!(scheme)
    
    #     factor1 = qFactor(scheme, node, edges[1])
    #     factor2 = qFactor(scheme, node, edges[2])
    #     factor3 = qFactor(scheme, node, edges[3])

    #     # Presetting marginals
    #     scheme.q_distributions[factor1] = GaussianDistribution()
    #     scheme.q_distributions[factor2] = GammaDistribution()
    #     scheme.q_distributions[factor3] = GaussianDistribution()
        
    #     # Univariate marginal
    #     calculateQDistribution!(node, factor1, scheme)
    #     @fact qDistribution(scheme, factor1) => GaussianDistribution(W=2.0, xi=0.0)
    #     # Univariate marginal
    #     calculateQDistribution!(node, factor2, scheme)
    #     @fact qDistribution(scheme, factor2) => GammaDistribution(a=1.0, b=2.0)
    #     # Univariate marginal
    #     calculateQDistribution!(node, factor3, scheme)
    #     @fact qDistribution(scheme, factor3) => GaussianDistribution(m=1.0, V=tiny())
    # end
    
    # context("Q distribution calculation for the structurally factorized GaussianNode") do
    #     (node, edges) = initializeGaussianNode(y_type=GaussianDistribution)
    #     scheme = InferenceScheme()
    #     factorize!(Set{Edge}([edges[3]]))
        
    #     factor1 = qFactor(scheme, node, edges[1])
    #     factor2 = qFactor(scheme, node, edges[2])
    #     @fact factor1 => factor2
    #     factor3 = qFactor(scheme, node, edges[3])

    #     # Presetting marginals
    #     scheme.q_distributions[factor1] = NormalGammaDistribution()
    #     scheme.q_distributions[factor3] = GaussianDistribution()
        
    #     # Joint marginal
    #     calculateQDistribution!(node, factor1, scheme)
    #     @fact qDistribution(scheme, factor1) => NormalGammaDistribution(m=0.0, beta=huge(), a=1.5, b=1.5)
    #     # Univariate marginal
    #     calculateQDistribution!(node, factor3, scheme)
    #     @fact qDistribution(scheme, factor3) => GaussianDistribution(W=2.0, xi=0.0)
    # end

        # # Mean field factorized Gaussian node
        # (node, edges) = initializeGaussianNode()
        # scheme = InferenceScheme()
        # factorize!(scheme)
        # setVagueQDistributions!(scheme)
        # sg_mean = subgraph(scheme, node.mean.edge)
        # sg_prec = subgraph(scheme, node.precision.edge)
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.mean, node.out)[1], qDistribution(scheme, node, sg_mean)) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.precision, node.out)[1], qDistribution(scheme, node, sg_prec)) => true

        # # Structurally factorized
        # (node, edges) = initializeGaussianNode()
        # scheme = InferenceScheme()
        # factorize!(node.out.edge)
        # setVagueQDistributions!(scheme)
        # sg_mean_prec = subgraph(scheme, node.mean.edge)
        # sg_out = subgraph(scheme, node.out.edge)
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.mean, node.out)[1], qDistribution(scheme, node, sg_mean_prec)) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.precision, node.out)[1], qDistribution(scheme, node, sg_mean_prec)) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.precision, node.mean)[1], node.precision.partner.message) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.mean, node.precision)[1], node.mean.partner.message) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.out, node.mean)[1], qDistribution(scheme, node, sg_out)) => true
        # @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.out, node.precision)[1], qDistribution(scheme, node, sg_out)) => true
