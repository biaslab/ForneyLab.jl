#####################
# Unit tests
#####################

facts("MvDelta unit tests") do
    context("MvDelta() should initialize a delta distribution") do
        @fact MvDelta().m --> [1.0]
        @fact typeof(MvDelta([2.0])) --> MvDelta{Float64, 1}
        @fact MvDelta([2.0]).m --> [2.0]
        @fact typeof(MvDelta([1.0, 2.0])) --> MvDelta{Float64, 2}
        @fact MvDelta([1.0, 2.0]).m --> [1.0, 2.0]
        @fact_throws MvDelta(reshape([2.0],1,1))
        @fact typeof(MvDelta()) --> MvDelta{Float64, 1}
        @fact pdf(MvDelta([2.0]), [2.0]) --> 1.0
        @fact pdf(MvDelta([2.0]), [2.1]) --> 0.0
    end

    context("MvDelta can be sampled") do
        @fact sample(MvDelta([2.0])) --> [2.0]
    end

    context("Product of two MvDeltas") do
        @fact MvDelta([2.0]) * MvDelta([2.0]) --> MvDelta([2.0])
        @fact_throws MvDelta([1.0]) * MvDelta([2.0])
        @fact_throws MethodError MvDelta([1.0]) * MvDelta([1.0; 1.0])
        @fact ForneyLab.prod!(MvDelta([2.0]), nothing) --> MvDelta([2.0])
    end

    context("Numbers and vectors should convert to MvDelta") do
        @fact convert(ProbabilityDistribution, [3.0, 4.0]) --> MvDelta([3.0, 4.0])
    end
end
