#####################
# Unit tests
#####################

facts("WishartDistribution unit tests") do
    context("WishartDistribution() should initialize a Wishart distribution") do
        dist = WishartDistribution(V=[2.0 1.0; 1.0 2.0], nu=3.0)
        @fact typeof(dist) --> WishartDistribution{2}
        @fact dist.V --> [2.0 1.0; 1.0 2.0]
        @fact dist.nu --> 3.0
        @fact mean(dist) --> 3.0*[2.0 1.0; 1.0 2.0]
        @fact var(dist) --> [24.0 15.0; 15.0 24.0]
        @fact var(WishartDistribution(V=diageye(3), nu=3.0)) --> [6.0 3.0 3.0; 3.0 6.0 3.0; 3.0 3.0 6.0]
    end

    context("vague() should initialize a vague (almost uninformative) Wishart distribution") do
        dist = vague(WishartDistribution{1})
        @fact dist.V --> [huge].'
        @fact dist.nu --> tiny

        dist = vague(WishartDistribution{2})
        @fact dist.V --> huge*eye(2)
        @fact dist.nu --> tiny
    end

    context("isProper() should indicate whether a Wishart distribution is proper") do
        @fact isProper(WishartDistribution()) --> true
        @fact isProper(WishartDistribution(V = [-1.0].', nu = 2.0)) --> false
        @fact isProper(WishartDistribution(V = [1.0].', nu = 0.0)) --> false
    end

    context("prod!() should yield correct result") do
        marg = WishartDistribution(V = [1.0].', nu=2.0) * WishartDistribution(V = [1.0].', nu=2.0)
        @fact marg.V --> roughly([0.5].')
        @fact marg.nu --> 2.0
    end
end

facts("WishartDistribution converts") do
    context("GammaDistribution should be convertible to WishartDistribution") do
        @fact convert(WishartDistribution, GammaDistribution(a=1.0, b=1.0)) --> WishartDistribution(V = [0.5].', nu = 2.0)
    end

    context("WishartDistribution should be convertible to GammaDistribution") do
        @fact convert(GammaDistribution, WishartDistribution(V = [0.5].', nu = 2.0)) --> GammaDistribution(a=1.0, b=1.0)
        @fact_throws convert(GammaDistribution, WishartDistribution(V=ones(2), nu=3.0))
    end
end

