facts("Partitioned unit tests") do
    context("Construction") do
        dd = Partitioned([Gaussian(), Gaussian()])
        @fact typeof(dd) --> Partitioned{Gaussian, 2}
        @fact dd.factors[1] --> Gaussian()
        @fact dd.factors[2] --> dd.factors[1]
        @fact typeof(Partitioned([vague(MvGaussian{2}), vague(MvGaussian{2})])) --> Partitioned{MvGaussian{2}, 2}
        @fact typeof(Partitioned([Gaussian()])) --> Partitioned{Gaussian, 1}
        @fact typeof(Partitioned([Gaussian(), Gamma()])) --> Partitioned{Union{Gaussian, Gamma}, 2}
        @fact typeof(Partitioned([vague(MvGaussian{2}), vague(MvGaussian{3})])) --> Partitioned{Union{MvGaussian{2}, MvGaussian{3}}, 2}
        dd = Partitioned([Gaussian(m=1.0, V=1.0), Gaussian()])
        @fact pdf(dd, [1.5;1.5]) --> roughly(pdf(Gaussian(m=1.0, V=1.0), 1.5) * pdf(Gaussian(), 1.5), atol=1e-6)
    end

    context("vague! and vague should be implemented") do
        dd = Partitioned([Gaussian(), Gaussian()])
        ForneyLab.vague!(dd)
        @fact dd.factors[1] --> vague(Gaussian)
        @fact dd.factors[2] --> vague(Gaussian)
        dtype = Partitioned{MvGaussian{2},3}
        vague_d = vague(dtype)
        @fact vague_d.factors[1] --> vague(MvGaussian{2})
        @fact vague_d.factors[2] --> vague_d.factors[1]
        d1 = Gaussian()
        d2 = Gamma()
        d = ForneyLab.vague!(Partitioned([d1, d2])) 
        @fact d --> Partitioned([vague(Gaussian), vague(Gamma)])
        @fact d1 --> vague(Gaussian)
        @fact d2 --> vague(Gamma)
        @fact_throws vague(Partitioned{Union{Gaussian, Gamma}})
    end

    context("mean should be implemented") do
        dd = Partitioned([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        @fact mean(dd) --> [1.0; 3.0]
        f1 = MvGaussian(m=2.0*ones(3),V=eye(3))
        f2 = MvGaussian(m=[6.;7.;8.],V=eye(3))
        @fact mean(Partitioned([f1;f2])) --> [2.;2.;2.;6.;7.;8.]
    end

    context("sample should be implemented") do
        s = sample(Partitioned([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)]))
        @fact typeof(s) --> Vector{Float64}
        @fact length(s) --> 2
    end

    context("== operator") do
        d1 = Partitioned([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        d2 = Partitioned([Gaussian(m=1.0,V=2.0), Gaussian(m=3.0,V=2.0)])
        d3 = Partitioned([Gaussian(m=1.0,V=2.0), Gaussian(m=3.1,V=2.0)])
        @fact d1 --> d2
        @fact (d1==d3) --> false
    end

    context("partitionedInboundTypesSize") do
        testfunc = ForneyLab.partitionedInboundTypesSize
        @fact testfunc([Gaussian; Message{Gamma}]) --> 0
        @fact testfunc([Partitioned{Gaussian,3}; Message{Gamma}]) --> 3
        @fact testfunc([Gaussian; Message{Partitioned{Gamma,2}}]) --> 2
    end

    context("stripPartitionedInboundTypes") do
        testfunc = ForneyLab.stripPartitionedInboundTypes
        @fact testfunc([Gaussian; Message{Gamma}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([Partitioned{Gaussian,3}; Message{Gamma}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([Gaussian; Message{Partitioned{Gamma,2}}]) --> [Gaussian; Message{Gamma}]
        @fact testfunc([Partitioned{Gaussian,2}; Message{Partitioned{Gamma,2}}]) --> [Gaussian; Message{Gamma}]
        @fact_throws testfunc([Partitioned{Gaussian,2}; Message{Partitioned{Gamma,3}}])
    end

    context("dimensions(::Partitioned)") do
        @fact dimensions(Partitioned([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])) --> 2
        @fact dimensions(Partitioned([Gaussian(m=1.0,V=1.0); Gamma(a=2.0,b=3.0)])) --> 2
        @fact dimensions(Partitioned([MvGaussian(m=zeros(2),V=eye(2)); MvGaussian(m=zeros(2),V=eye(2))])) --> 4
        @fact dimensions(Partitioned([MvGaussian(m=zeros(2),V=eye(2)); MvGaussian(m=zeros(3),V=eye(3))])) --> 5
    end

    context("prod() should yield correct result") do
        d1 = Partitioned([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])
        d2 = Partitioned([Gaussian(m=5.0,V=5.0); Gaussian(m=2.0,V=4.0)])
        marg = d1 * d2
        @fact typeof(marg) --> Partitioned{Gaussian, 2}
        @fact marg.factors[1] --> d1.factors[1] * d2.factors[1]
        @fact marg.factors[2] --> d1.factors[2] * d2.factors[2]
        @fact d1 * MvDelta(ones(2)) --> MvDelta(ones(2))
        d3 = Partitioned([MvGaussian(m=zeros(2),V=eye(2)); MvGaussian(m=zeros(2),V=eye(2))])
        @fact d3 * MvDelta(ones(4)) --> MvDelta(ones(4))
        @fact_throws DimensionMismatch d3 * MvDelta(ones(3))
        @fact_throws DomainError Partitioned([Gamma(); Gamma()]) * MvDelta(-1.*ones(2))
    end
end

facts("Partitioned integration tests") do
    context("Unrolling of partitioned inbounds") do
        FactorGraph()
        eq = EqualityNode()
        d1 = Partitioned([Gaussian(m=1.0,V=1.0); Gaussian(m=2.0,V=3.0)])
        d2 = Partitioned([Gaussian(m=5.0,V=5.0); Gaussian(m=2.0,V=4.0)])
        Edge(TerminalNode(d1), eq.interfaces[1])
        Edge(TerminalNode(d2), eq.interfaces[2])
        Edge(eq.interfaces[3], TerminalNode())

        algo = SumProduct(eq.interfaces[3])
        @fact algo.schedule[end].unrolling_factor --> 2
        execute(algo)
        out = eq.interfaces[3].message.payload

        @fact typeof(out) --> Partitioned{Gaussian, 2}
        @fact out.factors[1] --> d1.factors[1] * d2.factors[1]
        @fact out.factors[2] --> d1.factors[2] * d2.factors[2]
    end
end