facts("FactorizedDistribution unit tests") do
    context("Construction") do
        dd = FactorizedDistribution([GaussianDistribution(), GaussianDistribution()])
        @fact typeof(dd) --> FactorizedDistribution{GaussianDistribution, 2}
        @fact dd.factors[1] --> GaussianDistribution()
        @fact dd.factors[2] --> dd.factors[1]
        dd = FactorizedDistribution([vague(MvGaussianDistribution{2}), vague(MvGaussianDistribution{2})])
        @fact typeof(dd) --> FactorizedDistribution{MvGaussianDistribution{2}, 2}
        @fact_throws FactorizedDistribution([GaussianDistribution()])
        @fact_throws FactorizedDistribution([GaussianDistribution(), GammaDistribution()])
        @fact_throws FactorizedDistribution([vague(MvGaussianDistribution{2}), vague(MvGaussianDistribution{3})])
    end

    context("vague! and vague should be implemented") do
        dd = FactorizedDistribution([GaussianDistribution(), GaussianDistribution()])
        ForneyLab.vague!(dd)
        @fact dd.factors[1] --> vague(GaussianDistribution)
        @fact dd.factors[2] --> vague(GaussianDistribution)
        dtype = FactorizedDistribution{MvGaussianDistribution{2},3}
        vague_d = vague(dtype)
        @fact vague_d.factors[1] --> vague(MvGaussianDistribution{2})
        @fact vague_d.factors[2] --> vague_d.factors[1]
    end

    context("mean should be implemented") do
        dd = FactorizedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        @fact mean(dd) --> [1.0; 3.0]
        f1 = MvGaussianDistribution(m=2.0*ones(3),V=eye(3))
        f2 = MvGaussianDistribution(m=[6.;7.;8.],V=eye(3))
        @fact mean(FactorizedDistribution([f1;f2])) --> [2.;2.;2.;6.;7.;8.]
    end

    context("sample should be implemented") do
        s = sample(FactorizedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)]))
        @fact typeof(s) --> Vector{Float64}
        @fact length(s) --> 2
    end

    context("== operator") do
        d1 = FactorizedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        d2 = FactorizedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.0,V=2.0)])
        d3 = FactorizedDistribution([GaussianDistribution(m=1.0,V=2.0), GaussianDistribution(m=3.1,V=2.0)])
        @fact d1 --> d2
        @fact (d1==d3) --> false
    end
end
