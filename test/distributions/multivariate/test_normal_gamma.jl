#####################
# Unit tests
#####################

facts("NormalGammaDistribution unit tests") do
    context("NormalGammaDistribution() should initiatize a normal-gamma distribution") do
        dist = NormalGammaDistribution()
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> 1.0
        @fact dist.b --> 1.0
    end

    context("vague() should initialize an vague normal-gamma distribution") do
        dist = vague(NormalGammaDistribution)
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end

    context("prod!() should calculate product with a MvDeltaDistribution{Float64,2}") do
        @fact NormalGammaDistribution() * MvDeltaDistribution(ones(2)) --> MvDeltaDistribution(ones(2))
        @fact MvDeltaDistribution(ones(2)) * NormalGammaDistribution() --> MvDeltaDistribution(ones(2))
        @fact_throws MethodError MvDeltaDistribution(ones(3)) * NormalGammaDistribution()
        @fact_throws DomainError MvDeltaDistribution(-1.*ones(2)) * NormalGammaDistribution()
        @fact typeof(ForneyLab.prod!(NormalGammaDistribution(), MvDeltaDistribution(ones(2)), NormalGammaDistribution())) --> NormalGammaDistribution
        @fact mean(ForneyLab.prod!(MvDeltaDistribution(ones(2)), NormalGammaDistribution(), NormalGammaDistribution())) --> roughly([1.0;1.0])
        @fact var(ForneyLab.prod!(MvDeltaDistribution(ones(2)), NormalGammaDistribution(), NormalGammaDistribution())) --> roughly(zeros(2), atol=1e-8)
        @fact_throws MethodError ForneyLab.prod!(NormalGammaDistribution(), MvDeltaDistribution(ones(3)), NormalGammaDistribution())
        @fact_throws DomainError ForneyLab.prod!(NormalGammaDistribution(), MvDeltaDistribution(-1.*ones(2)), NormalGammaDistribution())
    end
end
