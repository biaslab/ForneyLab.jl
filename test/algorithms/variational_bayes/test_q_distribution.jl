# TODO: double check values

facts("Calculations for q distributions") do
    context("Q distribution calculation for naively factorized GaussianNode") do
        initializeGaussianNode()
        algo = VariationalBayes()
        prepare!(algo)
        f = algo.factorization
        qs = algo.q_distributions

        sm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        sp = f.edge_to_subgraph[n(:node).i[:precision].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        ForneyLab.calculateQDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), sm)].distribution --> GaussianDistribution(m=0.0, V=huge)
        ForneyLab.calculateQDistribution!(qs, n(:node), sp, f)
        @fact qs[(n(:node), sp)].distribution --> GammaDistribution(a=-0.999999999998, b=2.0e-12)
        ForneyLab.calculateQDistribution!(qs, n(:node), sm, f)
        @fact qs[(n(:node), sm)].distribution --> GaussianDistribution(m=0.0, V=5.0e11)
    end

    context("Q distribution calculation for the structurally factorized GaussianNode") do
        initializeGaussianNode()
        algo = VariationalBayes(Set{Edge}(Edge[ForneyLab.e(:edge3)]))
        prepare!(algo)
        f = algo.factorization
        qs = algo.q_distributions

        spm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        # Joint marginal
        ForneyLab.calculateQDistribution!(qs, n(:node), spm, f)
        @fact qs[(n(:node), spm)].distribution --> NormalGammaDistribution(m=0.0, beta=huge, a=0.500000000001, b=5.0e11)
        # Univariate marginal
        ForneyLab.calculateQDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), so)].distribution --> GaussianDistribution(xi=0.0, W=2e-12)
    end
end
