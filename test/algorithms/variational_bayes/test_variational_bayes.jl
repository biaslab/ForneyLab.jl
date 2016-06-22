facts("VariationalBayes should collect the proper inbound types as dependent on the recognition factorization") do
    # Mean field factorized Gaussian node
    initializeGaussianNode()

    rf = RecognitionFactorization()
    factorizeMeanField()
    initialize(eg(:edge1), vague(Gaussian))
    initialize(eg(:edge2), vague(Gamma))
    initialize(eg(:edge3), vague(Gaussian))

    algo = VariationalBayes()

    @fact algo.recognition_factorization.subgraphs[3].internal_schedule[2].inbound_types --> [Gaussian, Gamma, Void]
    @fact algo.recognition_factorization.subgraphs[2].internal_schedule[2].inbound_types --> [Gaussian, Void, Gaussian]
    @fact algo.recognition_factorization.subgraphs[1].internal_schedule[2].inbound_types --> [Void, Gamma, Gaussian]

    # Structurally factorized
    initializeGaussianNode()

    rf = RecognitionFactorization()
    factor([eg(:edge1), eg(:edge2)])
    factor(eg(:edge3))
    initialize([eg(:edge1), eg(:edge2)], vague(NormalGamma))
    initialize(eg(:edge3), vague(Gaussian))

    algo = VariationalBayes()
    
    @fact algo.recognition_factorization.subgraphs[2].internal_schedule[2].inbound_types --> [NormalGamma, NormalGamma, Void]
end
