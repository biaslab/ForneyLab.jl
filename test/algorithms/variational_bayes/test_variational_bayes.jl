facts("VariationalBayes collect inbound type tests") do
    context("VariationalBayes should collect the proper inbound types as dependent on the factorization") do
        # Mean field factorized Gaussian node
        initializeGaussianNode()

        algo = VariationalBayes(Dict(
            eg(:edge1) => GaussianDistribution,
            eg(:edge2) => GammaDistribution,
            eg(:edge3) => GaussianDistribution))

        @fact algo.factorization.factors[3].internal_schedule[2].inbound_types --> [GaussianDistribution, GammaDistribution, Void]
        @fact algo.factorization.factors[2].internal_schedule[2].inbound_types --> [GaussianDistribution, Void, GaussianDistribution]
        @fact algo.factorization.factors[1].internal_schedule[2].inbound_types --> [Void, GammaDistribution, GaussianDistribution]

        # Structurally factorized
        initializeGaussianNode()
        
        algo = VariationalBayes(Dict(
            eg(:edge*(1:2)).' => NormalGammaDistribution,
            eg(:edge3) => GaussianDistribution))
        
        @fact algo.factorization.factors[2].internal_schedule[2].inbound_types --> [NormalGammaDistribution, NormalGammaDistribution, Void]
    end
end