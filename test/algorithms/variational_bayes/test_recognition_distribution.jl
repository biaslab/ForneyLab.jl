# TODO: double check values

facts("Calculations for recognition distributions") do
    context("Recognition distribution calculation for naively factorized GaussianNode") do
        initializeGaussianNode()

        algo = VariationalBayes(Dict(
            eg(:edge1) => Gaussian,
            eg(:edge2) => Gamma,
            eg(:edge3) => Gaussian))
        
        prepare!(algo)
        f = algo.factorization
        qs = algo.recognition_distributions

        sm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        sp = f.edge_to_subgraph[n(:node).i[:precision].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        ForneyLab.calculateRecognitionDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), sm)].distribution --> Gaussian(m=0.0, V=huge)
        ForneyLab.calculateRecognitionDistribution!(qs, n(:node), sp, f)
        @fact qs[(n(:node), sp)].distribution --> Gamma(a=-0.999999999998, b=2.0e-12)
        ForneyLab.calculateRecognitionDistribution!(qs, n(:node), sm, f)
        @fact qs[(n(:node), sm)].distribution --> Gaussian(m=0.0, V=5.0e11)
    end

    context("Recognition distribution calculation for the structurally factorized GaussianNode") do
        initializeGaussianNode()

        algo = VariationalBayes(Dict(
            eg(:edge*(1:2)).' => NormalGamma,
            eg(:edge3) => Gaussian))
        
        prepare!(algo)
        f = algo.factorization
        qs = algo.recognition_distributions

        spm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        # Joint marginal
        ForneyLab.calculateRecognitionDistribution!(qs, n(:node), spm, f)
        @fact qs[(n(:node), spm)].distribution --> NormalGamma(m=0.0, beta=huge, a=0.500000000001, b=5.0e11)
        # Univariate marginal
        ForneyLab.calculateRecognitionDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), so)].distribution --> Gaussian(xi=0.0, W=2e-12)
    end
end
