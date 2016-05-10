facts("PartitionedDistribution unit tests") do
    context("Construction") do
        dd = PartitionedDistribution([Gaussian(), Gaussian()])
        @fact typeof(dd) --> PartitionedDistribution{Gaussian, 2}
        @fact dd.factors[1] --> Gaussian()
        @fact dd.factors[2] --> dd.factors[1]
        dd = PartitionedDistribution([vague(MvGaussian{2}), vague(MvGaussian{2})])
        @fact typeof(dd) --> PartitionedDistribution{MvGaussian{2}, 2}
        @fact_throws PartitionedDistribution([Gaussian()])
        @fact_throws PartitionedDistribution([Gaussian(), Gamma()])
        @fact_throws PartitionedDistribution([vague(MvGaussian{2}), vague(MvGaussian{3})])
    end

    context("vague! and vague should be implemented") do
        dd = PartitionedDistribution([Gaussian(), Gaussian()])
        ForneyLab.vague!(dd)
        @fact dd.factors[1] --> vague(Gaussian)
        @fact dd.factors[2] --> vague(Gaussian)
        dtype = PartitionedDistribution{MvGaussian{2},3}
        vague_d = vague(dtype)
        @fact vague_d.factors[1] --> vague(MvGaussian{2})
        @fact vague_d.factors[2] --> vague_d.factors[1]
    end

    context("mean should be implemented") do
        dd = PartitionedDistribution([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        @fact mean(dd) --> [1.0; 3.0]
        f1 = MvGaussian(m=2.0*ones(3),V=eye(3))
        f2 = MvGaussian(m=[6.;7.;8.],V=eye(3))
        @fact mean(PartitionedDistribution([f1;f2])) --> [2.;2.;2.;6.;7.;8.]
    end

    context("sample should be implemented") do
        s = sample(PartitionedDistribution([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)]))
        @fact typeof(s) --> Vector{Float64}
        @fact length(s) --> 2
    end

    context("== operator") do
        d1 = PartitionedDistribution([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        d2 = PartitionedDistribution([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        d3 = PartitionedDistribution([Gaussian(m=1.0,V=2.0), Gaussian(m=3.1,V=2.0)])
        @fact d1 --> d2
        @fact (d1==d3) --> false
    end

    context("partitionedInboundTypesSize") do
        testfunc = ForneyLab.partitionedInboundTypesSize
        @fact testfunc([Gaussian; Message{Gamma}]) --> 0
        @fact testfunc([PartitionedDistribution{Gaussian,3}; Message{Gamma}]) --> 3
        @fact testfunc([Gaussian; Message{PartitionedDistribution{Gamma,2}}]) --> 2
    end

    context("stripPartitionedInboundTypes") do
        testfunc = ForneyLab.stripPartitionedInboundTypes
        @fact testfunc([Gaussian; Message{Gamma}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([PartitionedDistribution{Gaussian,3}; Message{Gamma}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([Gaussian; Message{PartitionedDistribution{Gamma,2}}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([PartitionedDistribution{Gaussian,2}; Message{PartitionedDistribution{Gamma,2}}]) --> [Gaussian; Message{Gamma}]
        @fact_throws testfunc([PartitionedDistribution{Gaussian,2}; Message{PartitionedDistribution{Gamma,3}}])
    end

    context("dimensions(::PartitionedDistribution)") do
        @fact dimensions(PartitionedDistribution([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])) --> 2
        @fact dimensions(PartitionedDistribution([MvGaussian(m=zeros(2),V=eye(2)); MvGaussian(m=zeros(2),V=eye(2))])) --> 4
    end

    context("prod() should yield correct result") do
        d1 = PartitionedDistribution([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])
        d2 = PartitionedDistribution([Gaussian(m=5.0,V=5.0); Gaussian(m=2.0,V=4.0)])
        marg = d1 * d2
        @fact typeof(marg) --> PartitionedDistribution{Gaussian, 2}
        @fact marg.factors[1] --> d1.factors[1] * d2.factors[1]
        @fact marg.factors[2] --> d1.factors[2] * d2.factors[2]
        @fact d1 * MvDelta(ones(2)) --> MvDelta(ones(2))
        d3 = PartitionedDistribution([MvGaussian(m=zeros(2),V=eye(2)); MvGaussian(m=zeros(2),V=eye(2))])
        @fact d3 * MvDelta(ones(4)) --> MvDelta(ones(4))
        @fact_throws DimensionMismatch d3 * MvDelta(ones(3))
        @fact_throws DomainError PartitionedDistribution([Gamma(); Gamma()]) * MvDelta(-1.*ones(2))
    end
end

facts("PartitionedDistribution integration tests") do
    context("Unrolling of partitioned inbounds") do
        FactorGraph()
        eq = EqualityNode()
        d1 = PartitionedDistribution([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])
        d2 = PartitionedDistribution([Gaussian(m=5.0,V=5.0); Gaussian(m=2.0,V=4.0)])
        Edge(TerminalNode(d1), eq.interfaces[1])
        Edge(TerminalNode(d2), eq.interfaces[2])
        Edge(eq.interfaces[3], TerminalNode())

        algo = SumProduct(eq.interfaces[3])
        @fact algo.schedule[end].unrolling_factor --> 2
        execute(algo)
        out = eq.interfaces[3].message.payload

        @fact typeof(out) --> PartitionedDistribution{Gaussian, 2}
        @fact out.factors[1] --> d1.factors[1] * d2.factors[1]
        @fact out.factors[2] --> d1.factors[2] * d2.factors[2]
    end
end