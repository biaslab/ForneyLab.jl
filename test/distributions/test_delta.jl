#####################
# Unit tests
#####################

facts("DeltaDistribution unit tests") do
    context("DeltaDistribution() should initialize a delta distribution") do
        @fact DeltaDistribution().m => 1.0
        @fact DeltaDistribution(2.0).m => 2.0
        @fact DeltaDistribution([2.0]).m => [2.0]
        @fact DeltaDistribution(reshape([2.0],1,1)).m => reshape([2.0],1,1)
        @fact DeltaDistribution(:something).m => :something
    end

    context("DeltaDistribution can be sampled") do
        @fact sample(DeltaDistribution()) => DeltaDistribution()
    end

    context("There should be no such thing as vague(DeltaDistribution)") do
        @fact_throws vague(DeltaDistribution)
        @fact_throws vague(DeltaDistribution{Float64})
    end

    context("Numbers, symbols, and arrays should convert to DeltaDistribution") do
        @fact convert(ProbabilityDistribution, 3.0) => DeltaDistribution(3.0)
        @fact convert(ProbabilityDistribution, [3.0, 4.0]) => DeltaDistribution([3.0, 4.0])
        @fact convert(ProbabilityDistribution, :something) => DeltaDistribution(:something)
    end
end