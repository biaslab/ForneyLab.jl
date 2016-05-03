facts("PartitionedDistribution unit tests") do
    context("Construction") do
        dd = PartitionedDistribution([GaussianDistribution(), GaussianDistribution()])
        @fact typeof(dd) --> PartitionedDistribution{GaussianDistribution, 2}
        @fact dd.factors[1] --> GaussianDistribution()
        @fact dd.factors[2] --> dd.factors[1]
        dd = PartitionedDistribution([vague(MvGaussianDistribution{2}), vague(MvGaussianDistribution{2})])
        @fact typeof(dd) --> PartitionedDistribution{MvGaussianDistribution{2}, 2}
        @fact_throws PartitionedDistribution([GaussianDistribution()])
        @fact_throws PartitionedDistribution([GaussianDistribution(), GammaDistribution()])
        @fact_throws PartitionedDistribution([vague(MvGaussianDistribution{2}), vague(MvGaussianDistribution{3})])
    end

    context("vague! and vague should be implemented") do
        dd = PartitionedDistribution([GaussianDistribution(), GaussianDistribution()])
        ForneyLab.vague!(dd)
        @fact dd.factors[1] --> vague(GaussianDistribution)
        @fact dd.factors[2] --> vague(GaussianDistribution)
        dtype = PartitionedDistribution{MvGaussianDistribution{2},3}
        vague_d = vague(dtype)
        @fact vague_d.factors[1] --> vague(MvGaussianDistribution{2})
        @fact vague_d.factors[2] --> vague_d.factors[1]
    end

    context("mean should be implemented") do
        dd = PartitionedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        @fact mean(dd) --> [1.0; 3.0]
        f1 = MvGaussianDistribution(m=2.0*ones(3),V=eye(3))
        f2 = MvGaussianDistribution(m=[6.;7.;8.],V=eye(3))
        @fact mean(PartitionedDistribution([f1;f2])) --> [2.;2.;2.;6.;7.;8.]
    end

    context("sample should be implemented") do
        s = sample(PartitionedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)]))
        @fact typeof(s) --> Vector{Float64}
        @fact length(s) --> 2
    end

    context("== operator") do
        d1 = PartitionedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        d2 = PartitionedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        d3 = PartitionedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.1,V=2.0)])
        @fact d1 --> d2
        @fact (d1==d3) --> false
    end

    context("partitionedInboundTypesSize") do
        testfunc = ForneyLab.partitionedInboundTypesSize
        @fact testfunc([GaussianDistribution; Message{GammaDistribution}]) --> 0
        @fact testfunc([PartitionedDistribution{GaussianDistribution,3}; Message{GammaDistribution}]) --> 3
        @fact testfunc([GaussianDistribution; Message{PartitionedDistribution{GammaDistribution,2}}]) --> 2
    end

    context("stripPartitionedInboundTypes") do
        testfunc = ForneyLab.stripPartitionedInboundTypes
        @fact testfunc([GaussianDistribution; Message{GammaDistribution}]) --> [GaussianDistribution; Message{GammaDistribution}]
        @fact testfunc([PartitionedDistribution{GaussianDistribution,3}; Message{GammaDistribution}]) --> [GaussianDistribution; Message{GammaDistribution}]
        @fact testfunc([GaussianDistribution; Message{PartitionedDistribution{GammaDistribution,2}}]) --> [GaussianDistribution; Message{GammaDistribution}]
        @fact testfunc([PartitionedDistribution{GaussianDistribution,2}; Message{PartitionedDistribution{GammaDistribution,2}}]) --> [GaussianDistribution; Message{GammaDistribution}]
        @fact_throws testfunc([PartitionedDistribution{GaussianDistribution,2}; Message{PartitionedDistribution{GammaDistribution,3}}])
    end

    context("calculateMarginal") do
        d1 = PartitionedDistribution([GaussianDistribution(m=1.0,V=1.0); GaussianDistribution(m=2.0,V=3.0)])
        d2 = PartitionedDistribution([GaussianDistribution(m=5.0,V=5.0); GaussianDistribution(m=2.0,V=4.0)])
        marg = calculateMarginal(d1, d2)
        @fact typeof(marg) --> PartitionedDistribution{GaussianDistribution, 2}
        @fact marg.factors[1] --> calculateMarginal(d1.factors[1], d2.factors[1])
        @fact marg.factors[2] --> calculateMarginal(d1.factors[2], d2.factors[2])
    end
end

facts("PartitionedDistribution integration tests") do
    context("Unrolling of partitioned inbounds") do
        FactorGraph()
        eq = EqualityNode()
        d1 = PartitionedDistribution([GaussianDistribution(m=1.0,V=1.0); GaussianDistribution(m=2.0,V=3.0)])
        d2 = PartitionedDistribution([GaussianDistribution(m=5.0,V=5.0); GaussianDistribution(m=2.0,V=4.0)])
        Edge(TerminalNode(d1), eq.interfaces[1])
        Edge(TerminalNode(d2), eq.interfaces[2])
        Edge(eq.interfaces[3], TerminalNode())

        algo = SumProduct(eq.interfaces[3])
        @fact algo.schedule[end].unrolling_factor --> 2
        execute(algo)
        out = eq.interfaces[3].message.payload

        @fact typeof(out) --> PartitionedDistribution{GaussianDistribution, 2}
        @fact out.factors[1] --> calculateMarginal(d1.factors[1], d2.factors[1])
        @fact out.factors[2] --> calculateMarginal(d1.factors[2], d2.factors[2])
    end
end