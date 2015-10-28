#####################
# Unit tests
#####################

facts("MvDeltaDistribution unit tests") do
    context("MvDeltaDistribution() should initialize a delta distribution") do
        @fact MvDeltaDistribution().m => [1.0]
        @fact MvDeltaDistribution([2.0]).m => [2.0]
        @fact MvDeltaDistribution([1.0, 2.0]).m => [1.0, 2.0]
        @fact_throws MvDeltaDistribution(reshape([2.0],1,1))
        @fact typeof(MvDeltaDistribution()) => MvDeltaDistribution{Float64}
    end

    context("MvDeltaDistribution can be sampled") do
        @fact sample(MvDeltaDistribution([2.0])) => [2.0]
    end

    context("There should be no such thing as vague(MvDeltaDistribution)") do
        @fact_throws vague(MvDeltaDistribution)
        @fact_throws vague(MvDeltaDistribution{Float64})
    end

    context("Marginal calculation for two MvDeltaDistributions") do
        @fact calculateMarginal(MvDeltaDistribution([2.0]), MvDeltaDistribution([2.0])) => MvDeltaDistribution([2.0])
        @fact calculateMarginal(MvDeltaDistribution([1.0]), MvDeltaDistribution([2.0])) => MvDeltaDistribution([0.0])
    end

    context("Numbers and vectors should convert to MvDeltaDistribution") do
        @fact convert(ProbabilityDistribution, [3.0, 4.0]) => MvDeltaDistribution([3.0, 4.0])
    end
end