#####################
# Unit tests
#####################

facts("BernoulliDistribution unit tests") do
    context("construction") do
        @fact BernoulliDistribution().p --> 0.5
        @fact BernoulliDistribution(0.8).p --> 0.8
    end

    context("required methods") do
        @fact vague(BernoulliDistribution) --> BernoulliDistribution(0.5)
        @fact mean(BernoulliDistribution(0.7)) --> 0.7
        @fact var(BernoulliDistribution(0.7)) --> roughly(0.7*0.3)
        @fact typeof(sample(BernoulliDistribution())) --> Bool
        @fact sample(BernoulliDistribution(0.0)) --> false
        @fact sample(BernoulliDistribution(1.0)) --> true
        @fact (BernoulliDistribution() == BernoulliDistribution()) --> true
        @fact (BernoulliDistribution(0.3) == BernoulliDistribution(0.7)) --> false
    end

    context("prod!") do
        @fact (BernoulliDistribution(0.2) * BernoulliDistribution(0.4)).p --> roughly(1/7)
        @fact_throws BernoulliDistribution(0.0) * BernoulliDistribution(1.0)
        @fact BernoulliDistribution(0.2) * DeltaDistribution(true) --> DeltaDistribution(true)
        @fact DeltaDistribution(false) * BernoulliDistribution(0.2) --> DeltaDistribution(false)
        @fact_throws BernoulliDistribution(0.0) * DeltaDistribution(true)
        @fact_throws BernoulliDistribution(0.2) * DeltaDistribution(3.0)
        @fact ForneyLab.prod!(BernoulliDistribution(0.2), DeltaDistribution(true), BernoulliDistribution()) --> BernoulliDistribution(1.0)
        @fact ForneyLab.prod!(DeltaDistribution(true), BernoulliDistribution(0.2), BernoulliDistribution()) --> BernoulliDistribution(1.0)
        @fact_throws ForneyLab.prod!(DeltaDistribution(false), BernoulliDistribution(1.0), BernoulliDistribution())
    end
end

#####################
# Integration tests
#####################

facts("BernoulliDistribution integration tests") do
    context("Convert DeltaDistribution to BernoulliDistribution") do
        @fact convert(BernoulliDistribution, DeltaDistribution(false)) --> BernoulliDistribution(0.0)
        @fact convert(BernoulliDistribution, DeltaDistribution(true)) --> BernoulliDistribution(1.0)
        @fact_throws convert(BernoulliDistribution, DeltaDistribution(2.0))
    end
end
