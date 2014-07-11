#####################
# Unit tests
#####################

facts("GammaDistribution unit tests") do
    context("GammaDistribution() should initiatize a gamma distribution") do
        dist = GammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
    end

    context("uninformative() should initialize an uninformative gamma distribution") do
        dist = uninformative(GammaDistribution)
        @fact dist.a => 0.999
        @fact dist.b => 0.001
    end
end