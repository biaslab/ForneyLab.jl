#####################
# Unit tests
#####################

facts("NormalGammaDistribution unit tests") do
    context("NormalGammaDistribution() should initiatize a normal-gamma distribution") do
        dist = NormalGammaDistribution()
        @fact dist.m => 0.0
        @fact dist.W => 1.0
        @fact dist.a => 1.0
        @fact dist.b => 1.0
    end

    context("uninformative() should initialize an uninformative normal-gamma distribution") do
        dist = uninformative(NormalGammaDistribution)
        @fact dist.m => 0.0
        @fact dist.W => 1.0
        @fact dist.a => 0.999
        @fact dist.b => 0.001
    end
end
