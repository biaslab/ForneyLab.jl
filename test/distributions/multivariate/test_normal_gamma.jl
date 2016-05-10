#####################
# Unit tests
#####################

facts("NormalGamma unit tests") do
    context("NormalGamma() should initiatize a normal-gamma distribution") do
        dist = NormalGamma()
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> 1.0
        @fact dist.b --> 1.0
    end

    context("vague() should initialize an vague normal-gamma distribution") do
        dist = vague(NormalGamma)
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> tiny
        @fact dist.b --> tiny
    end

    context("prod!() should calculate product with a MvDelta{Float64,2}") do
        @fact NormalGamma() * MvDelta(ones(2)) --> MvDelta(ones(2))
        @fact MvDelta(ones(2)) * NormalGamma() --> MvDelta(ones(2))
        @fact_throws MethodError MvDelta(ones(3)) * NormalGamma()
        @fact_throws DomainError MvDelta(-1.*ones(2)) * NormalGamma()
        @fact typeof(ForneyLab.prod!(NormalGamma(), MvDelta(ones(2)), NormalGamma())) --> NormalGamma
        @fact mean(ForneyLab.prod!(MvDelta(ones(2)), NormalGamma(), NormalGamma())) --> roughly([1.0;1.0])
        @fact var(ForneyLab.prod!(MvDelta(ones(2)), NormalGamma(), NormalGamma())) --> roughly(zeros(2), atol=1e-8)
        @fact_throws MethodError ForneyLab.prod!(NormalGamma(), MvDelta(ones(3)), NormalGamma())
        @fact_throws DomainError ForneyLab.prod!(NormalGamma(), MvDelta(-1.*ones(2)), NormalGamma())
    end
end
