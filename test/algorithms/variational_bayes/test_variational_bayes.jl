facts("VariationalBayes should collect the proper inbound types as dependent on the recognition factorization") do
    # Mean field factorized Gaussian node
    initializeGaussianNode()

    rf = RecognitionFactorization()
    factorizeMeanField()
    setRecognitionDistribution(eg(:edge1), Gaussian)
    setRecognitionDistribution(eg(:edge2), Gamma)
    setRecognitionDistribution(eg(:edge3), Gaussian)

    algo = VariationalBayes()

    @fact algo.recognition_factorization.subgraphs[3].internal_schedule[2].inbound_types --> [Gaussian, Gamma, Void]
    @fact algo.recognition_factorization.subgraphs[2].internal_schedule[2].inbound_types --> [Gaussian, Void, Gaussian]
    @fact algo.recognition_factorization.subgraphs[1].internal_schedule[2].inbound_types --> [Void, Gamma, Gaussian]

    # Structurally factorized
    initializeGaussianNode()

    rf = RecognitionFactorization()
    addFactor([eg(:edge1), eg(:edge2)])
    addFactor(eg(:edge3))
    setRecognitionDistribution([eg(:edge1), eg(:edge2)], NormalGamma)
    setRecognitionDistribution(eg(:edge3), Gaussian)

    algo = VariationalBayes()
    
    @fact algo.recognition_factorization.subgraphs[2].internal_schedule[2].inbound_types --> [NormalGamma, NormalGamma, Void]
end
