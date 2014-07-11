#####################
# Unit tests
#####################

facts("InverseGammaDistribution unit tests") do
    context("InverseGammaDistribution() should initiatize an inverse gamma distribution") do
        dist = InverseGammaDistribution(a=2.0, b=0.5)
        @fact dist.a => 2.0
        @fact dist.b => 0.5
    end

    context("uninformative() should initialize an uninformative inverse gamma distribution") do
        dist = uninformative(InverseGammaDistribution)
        @fact dist.a => -0.999
        @fact dist.b => 0.001
    end
end