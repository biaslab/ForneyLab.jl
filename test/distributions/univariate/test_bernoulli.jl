#####################
# Unit tests
#####################

facts("Bernoulli unit tests") do
    context("construction") do
        @fact Bernoulli().p --> 0.5
        @fact Bernoulli(0.8).p --> 0.8
    end

    context("required methods") do
        @fact vague(Bernoulli) --> Bernoulli(0.5)
        @fact mean(Bernoulli(0.7)) --> 0.7
        @fact var(Bernoulli(0.7)) --> roughly(0.7*0.3)
        @fact pdf(Bernoulli(0.6), true) --> 0.6
        @fact pdf(Bernoulli(0.6), false) --> 0.4
        @fact typeof(sample(Bernoulli())) --> Bool
        @fact sample(Bernoulli(0.0)) --> false
        @fact sample(Bernoulli(1.0)) --> true
        @fact (Bernoulli() == Bernoulli()) --> true
        @fact (Bernoulli(0.3) == Bernoulli(0.7)) --> false
    end

    context("prod!") do
        @fact (Bernoulli(0.2) * Bernoulli(0.4)).p --> roughly(1/7)
        @fact_throws Bernoulli(0.0) * Bernoulli(1.0)
        @fact Bernoulli(0.2) * Delta(true) --> Delta(true)
        @fact Delta(false) * Bernoulli(0.2) --> Delta(false)
        @fact_throws Bernoulli(0.0) * Delta(true)
        @fact_throws Bernoulli(0.2) * Delta(3.0)
        @fact ForneyLab.prod!(Bernoulli(0.2), Delta(true), Bernoulli()) --> Bernoulli(1.0)
        @fact ForneyLab.prod!(Delta(true), Bernoulli(0.2), Bernoulli()) --> Bernoulli(1.0)
        @fact_throws ForneyLab.prod!(Delta(false), Bernoulli(1.0), Bernoulli())
    end
end

#####################
# Integration tests
#####################

facts("Bernoulli integration tests") do
    context("Convert Delta to Bernoulli") do
        @fact convert(Bernoulli, Delta(false)) --> Bernoulli(0.0)
        @fact convert(Bernoulli, Delta(true)) --> Bernoulli(1.0)
        @fact_throws convert(Bernoulli, Delta(2.0))
    end
end
