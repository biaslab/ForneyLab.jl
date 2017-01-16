#####################
# Unit tests
#####################

facts("Delta unit tests") do
    context("Delta() should initialize a delta distribution") do
        @fact Delta().m --> 1.0
        @fact Delta(2.0).m --> 2.0
        @fact typeof(Delta()) --> Delta{Float64}
        @fact pdf(Delta(2.0), 2.0) --> 1.0
        @fact pdf(Delta(2.0), 2.1) --> 0.0
    end

    context("Delta can be sampled") do
        @fact sample(Delta(2.0)) --> 2.0
    end

    context("prod! involving Deltas") do
        @fact Delta(2.0) * Delta(2.0) --> Delta(2.0)
        @fact_throws Delta(1.0) * Delta(2.0)
        @fact ForneyLab.prod!(Delta(2.0), nothing) --> Delta(2.0)
    end


    context("unsafeLogMean() should return correct result") do
        @fact ForneyLab.unsafeLogMean(Delta(2.0)) --> log(2)
    end

    context("Numbers should convert to Delta") do
        @fact convert(ProbabilityDistribution, 3.0) --> Delta(3.0)
    end
end
