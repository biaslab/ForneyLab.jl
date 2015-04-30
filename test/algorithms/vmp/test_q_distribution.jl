facts("Calculations for q distributions") do
    context("Q distribution calculation for naively factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=Float64)
        algo = VMP.Algorithm()
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    

        sm = f.edge_to_subgraph[node.mean.edge]
        sp = f.edge_to_subgraph[node.precision.edge]
        so = f.edge_to_subgraph[node.out.edge]

        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, so, f)
        @fact qs[(node, sm)].distribution.m[1] => 0.0
        @fact qs[(node, sm)].distribution.V[1,1] => huge()

        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, sp, f)
        @fact qs[(node, sp)].distribution => GammaDistribution(a=1.0, b=2.0)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, sm, f)
        @fact qs[(node, so)].distribution.m[1] => 1.0
        @fact qs[(node, so)].distribution.V[1,1] => tiny()
    end
    
    context("Q distribution calculation for the structurally factorized GaussianNode") do
        (node, edges) = initializeGaussianNode(y_type=GaussianDistribution)
        algo = VMP.Algorithm(Set{Edge}({edges[3]}))
        f = algo.fields[:factorization]
        qs = algo.fields[:q_distributions]    
        
        spm = f.edge_to_subgraph[node.mean.edge]
        so = f.edge_to_subgraph[node.out.edge]

        # Joint marginal
        VMP.calculateQDistribution!(qs, node, spm, f)
        @fact qs[(node, spm)].distribution => NormalGammaDistribution(m=0.0, beta=huge(), a=1.5, b=5.00000000001e11)
        # Univariate marginal
        VMP.calculateQDistribution!(qs, node, so, f)
        @fact qs[(node, so)].distribution => GaussianDistribution(W=2.0, xi=0.0)
    end
end