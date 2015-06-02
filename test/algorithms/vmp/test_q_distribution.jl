facts("Calculations for q distributions") do
    context("Q distribution calculation for naively factorized GaussianNode") do
        initializeGaussianNode(y_type=Float64)
        algo = VMP.Algorithm()
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    

        sm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        sp = f.edge_to_subgraph[n(:node).i[:precision].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        # Univariate marginal
        VMP.calculateQDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), sm)].distribution.m[1] => 0.0
        @fact qs[(n(:node), sm)].distribution.V[1,1] => huge()

        # Univariate marginal
        VMP.calculateQDistribution!(qs, n(:node), sp, f)
        @fact qs[(n(:node), sp)].distribution => GammaDistribution(a=1.0, b=2.0)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, n(:node), sm, f)
        @fact qs[(n(:node), so)].distribution.m[1] => 1.0
        @fact qs[(n(:node), so)].distribution.V[1,1] => tiny()
    end
    
    context("Q distribution calculation for the structurally factorized GaussianNode") do
        initializeGaussianNode(y_type=GaussianDistribution)
        algo = VMP.Algorithm(Set{Edge}({e(:edge3)}))
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    
        
        spm = f.edge_to_subgraph[n(:node).i[:mean].edge]
        so = f.edge_to_subgraph[n(:node).i[:out].edge]

        # Joint marginal
        VMP.calculateQDistribution!(qs, n(:node), spm, f)
        @fact qs[(n(:node), spm)].distribution => NormalGammaDistribution(m=0.0, beta=huge(), a=1.5, b=5.00000000001e11)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, n(:node), so, f)
        @fact qs[(n(:node), so)].distribution => GaussianDistribution(W=2.0, xi=0.0)
    end
end