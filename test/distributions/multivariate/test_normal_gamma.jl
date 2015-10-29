#####################
# Unit tests
#####################

facts("NormalGammaDistribution unit tests") do
    context("NormalGammaDistribution() should initiatize a normal-gamma distribution") do
        dist = NormalGammaDistribution()
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> 1.0
        @fact dist.b --> 1.0
    end

    context("vague() should initialize an vague normal-gamma distribution") do
        dist = vague(NormalGammaDistribution)
        @fact dist.m --> 0.0
        @fact dist.beta --> 1.0
        @fact dist.a --> 1.0-tiny
        @fact dist.b --> tiny
    end
end
