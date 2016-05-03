#####################
# Unit tests
#####################

facts("StudentsTDistribution unit tests") do
    context("StudentsTDistribution() should initiatize a student's t-distribution") do
        dist = StudentsTDistribution()
        @fact dist.m --> 0.0
        @fact dist.lambda --> 1.0
        @fact dist.nu --> huge
        @fact mean(StudentsTDistribution(m=0.0, lambda=1.0, nu=3.0)) --> 0.0
        @fact var(StudentsTDistribution(m=0.0, lambda=1.0, nu=3.0)) --> 3.0
        @fact isnan(mean(StudentsTDistribution(m=0.0, lambda=1.0, nu=1.0))) --> true
        @fact isnan(var(StudentsTDistribution(m=0.0, lambda=1.0, nu=0.5))) --> true
        @fact vague(StudentsTDistribution) --> StudentsTDistribution(m=0.0, lambda=tiny, nu=4.0)
    end

    context("prod!() should yield correct result") do
        @fact DeltaDistribution(0.5) * StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0) --> DeltaDistribution(0.5)
        @fact typeof(ForneyLab.prod!(DeltaDistribution(0.5), StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0), StudentsTDistribution())) --> StudentsTDistribution
        @fact mean(ForneyLab.prod!(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0), DeltaDistribution(0.5), StudentsTDistribution())) --> roughly(0.5)
        @fact var(ForneyLab.prod!(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0), DeltaDistribution(0.5), StudentsTDistribution())) --> less_than(1e-6)
    end

    context("Product of Gaussian and student's t-distribution should yield moment-matched Gaussian") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(GaussianDistribution()), Message(StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0)), nothing],
                                GaussianDistribution(m=0.5, W=2.0),
                                ForneyLab.sumProductRule!,
                                MomentMatching)
    end
end
