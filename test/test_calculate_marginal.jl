facts("Marginal calculation integration tests") do
    context("Marginal calculation for two GammaDistributions") do
        @fact calculateMarginal(GammaDistribution(a=1.0, b=2.0), GammaDistribution(a=3.0, b=4.0)) => GammaDistribution(a=3.0, b=6.0)
    end

    context("Marginal calculation for two InverseGammaDistributions") do
        @fact calculateMarginal(InverseGammaDistribution(a=1.0, b=2.0), InverseGammaDistribution(a=3.0, b=4.0)) => InverseGammaDistribution(a=5.0, b=6.0)
    end

    context("Marginal calculation for two GaussianDistributions") do
        @fact calculateMarginal(GaussianDistribution(xi=1.0, W=2.0), GaussianDistribution(xi=3.0, W=4.0)) => GaussianDistribution(xi=4.0, W=6.0)
    end

    context("Marginal calculation for the combination of a Gaussian and student's t-distribution") do
        @fact calculateMarginal(GaussianDistribution(m=0.0, W=1.0), StudentsTDistribution(m=0.0, W=1.0, nu=1.0)) => GaussianDistribution(m=0.0, W=3.0)
        @fact calculateMarginal(StudentsTDistribution(m=0.0, W=1.0, nu=1.0), GaussianDistribution(m=0.0, W=1.0)) => GaussianDistribution(m=0.0, W=3.0)
    end

    context("Q distribution calculation for naively factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=Float64)
        scheme = InferenceScheme()
        factorize!(scheme)
    
        factor1 = qFactor(scheme, node, edges[1])
        factor2 = qFactor(scheme, node, edges[2])
        factor3 = qFactor(scheme, node, edges[3])

        # Presetting marginals
        scheme.q_distributions[factor1] = GaussianDistribution()
        scheme.q_distributions[factor2] = GammaDistribution()
        scheme.q_distributions[factor3] = GaussianDistribution()
        
        # Univariate marginal
        calculateQDistribution!(node, factor1, scheme)
        @fact qDistribution(factor1) => GaussianDistribution(W=2.0, xi=0.0)
        # Univariate marginal
        calculateQDistribution!(node, factor2, scheme)
        @fact qDistribution(factor2) => GammaDistribution(a=1.0, b=2.0)
        # Univariate marginal
        calculateQDistribution!(node, factor3, scheme)
        @fact qDistribution(factor3) => GaussianDistribution(m=1.0, V=tiny())
    end
    
    context("Q distribution calculation for the structurally factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=GaussianDistribution)
        scheme = InferenceScheme()
        factorize!(Set{Edge}([edges[3]]))
        
        factor1 = qFactor(scheme, node, edges[1])
        factor2 = qFactor(scheme, node, edges[2])
        @fact factor1 => factor2
        factor3 = qFactor(scheme, node, edges[3])

        # Presetting marginals
        scheme.q_distributions[factor1] = NormalGammaDistribution()
        scheme.q_distributions[factor3] = GaussianDistribution()
        
        # Joint marginal
        calculateQDistribution!(node, factor1, scheme)
        @fact qDistribution(factor1) => NormalGammaDistribution(m=0.0, beta=huge(), a=1.5, b=1.5)
        # Univariate marginal
        calculateQDistribution!(node, factor3, scheme)
        @fact qDistribution(factor3) => GaussianDistribution(W=2.0, xi=0.0)
    end
end