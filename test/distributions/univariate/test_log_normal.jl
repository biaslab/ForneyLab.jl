#####################
# Unit tests
#####################

facts("LogNormalDistribution unit tests") do
    context("LogNormalDistribution() should initiatize a log-normal distribution") do
        dist = LogNormalDistribution(m=2.0, s=0.5)
        @fact dist.m --> 2.0
        @fact dist.s --> 0.5
        @fact mean(dist) --> exp(dist.m + 0.5*dist.s)
        @fact var(dist) --> (exp(dist.s) - 1)*exp(2*dist.m + dist.s)
        @fact isnan(mean(LogNormalDistribution(m=0.0, s=-1.0))) --> true
    end

    context("vague() should initialize a vague (almost uninformative) log-normal distribution") do
        dist = vague(LogNormalDistribution)
        @fact dist.m --> 0.0
        @fact dist.s --> huge
    end

    context("approximateWithGamma() should approximate a log-normal distribution with moment matching") do
        d = approximateWithLogNormal(GammaDistribution(a=1.0, b=5.0))
        @fact d.m --> -1.956011502714073
        @fact d.s --> 0.6931471805599453
    end

    context("approximateWithLogNormal() should approximate a gamma distribution with moment matching") do
        d = approximateWithGamma(LogNormalDistribution(m=-1.956011502714073, s=0.6931471805599453))
        @fact d.a --> 1.0
        @fact d.b --> 5.0
    end
end