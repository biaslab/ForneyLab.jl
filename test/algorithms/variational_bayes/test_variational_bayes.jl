facts("VariationalBayes collect inbound type tests") do
    context("VariationalBayes should collect the proper inbound types as dependent on the factorization") do
        # Mean field factorized Gaussian node
        initializeGaussianNode()

        algo = VariationalBayes(Dict(
            eg(:edge1) => Gaussian,
            eg(:edge2) => Gamma,
            eg(:edge3) => Gaussian))

        @fact algo.factorization.factors[3].internal_schedule[2].inbound_types --> [Gaussian, Gamma, Void]
        @fact algo.factorization.factors[2].internal_schedule[2].inbound_types --> [Gaussian, Void, Gaussian]
        @fact algo.factorization.factors[1].internal_schedule[2].inbound_types --> [Void, Gamma, Gaussian]

        # Structurally factorized
        initializeGaussianNode()
        
        algo = VariationalBayes(Dict(
            eg(:edge*(1:2)).' => NormalGamma,
            eg(:edge3) => Gaussian))
        
        @fact algo.factorization.factors[2].internal_schedule[2].inbound_types --> [NormalGamma, NormalGamma, Void]
    end
end