facts("Calculations for q distributions") do
    context("Q distribution calculation for naively factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=Float64)
        algo = VMP.Algorithm()
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    

        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, f.factors[1], f)
        @fact qs[(node, f.factors[1])].distribution => GaussianDistribution(m=1.0, V=tiny())
        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, f.factors[2], f)
        @fact qs[(node, f.factors[2])].distribution => GammaDistribution(a=1.0, b=2.0)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, f.factors[3], f)
        @fact qs[(node, f.factors[3])].distribution => GaussianDistribution(W=2.0, xi=0.0)
    end
    
    context("Q distribution calculation for the structurally factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=GaussianDistribution)
        algo = VMP.Algorithm([Set{Edge}([edges[3]])])
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    
        
        # Joint marginal
        VMP.calculateQDistribution!(qs, node, f.factors[1], f)
        @fact qs[(node, f.factors[1])].distribution => NormalGammaDistribution(m=0.0, beta=huge(), a=1.5, b=5.00000000001e11)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, f.factors[2], f)
        @fact qs[(node, f.factors[2])].distribution => GaussianDistribution(W=2.0, xi=0.0)
    end
end