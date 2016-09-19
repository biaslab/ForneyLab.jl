#####################
# Unit tests
#####################

facts("Wishart unit tests") do
    context("Wishart() should initialize a Wishart distribution") do
        dist = Wishart(V=[2.0 1.0; 1.0 2.0], nu=3.0)
        @fact typeof(dist) --> Wishart{2}
        @fact dist.V --> [2.0 1.0; 1.0 2.0]
        @fact dist.nu --> 3.0
        @fact mean(dist) --> 3.0*[2.0 1.0; 1.0 2.0]
        @fact var(dist) --> [24.0 15.0; 15.0 24.0]
        @fact var(Wishart(V=diageye(3), nu=3.0)) --> [6.0 3.0 3.0; 3.0 6.0 3.0; 3.0 3.0 6.0]
        @fact dimensions(Wishart{2}) --> (2, 2)
        @fact dimensions(Wishart(V=[2.0 1.0; 1.0 2.0], nu=3.0)) --> (2, 2)
    end

    context("vague() should initialize a vague (almost uninformative) Wishart distribution") do
        dist = vague(Wishart{1})
        @fact dist.V --> [huge].'
        @fact dist.nu --> tiny

        dist = vague(Wishart{2})
        @fact dist.V --> huge*eye(2)
        @fact dist.nu --> tiny
    end

    context("isProper() should indicate whether a Wishart distribution is proper") do
        @fact isProper(Wishart()) --> true
        @fact isProper(Wishart(V = [-1.0].', nu = 2.0)) --> false
        @fact isProper(Wishart(V = [1.0].', nu = 0.0)) --> false
    end

    context("prod!() should yield correct result") do
        marg = Wishart(V = [1.0].', nu=2.0) * Wishart(V = [1.0].', nu=2.0)
        @fact marg.V --> roughly([0.5].')
        @fact marg.nu --> 2.0
    end

    context("unsafeDetLogMean() should return correct result") do
        @fact ForneyLab.unsafeDetLogMean(Wishart(V=[1.0].', nu=1.0)) --> digamma(0.5) + log(2)
        @fact ForneyLab.unsafeDetLogMean(Wishart(V=eye(2), nu=2.0)) --> digamma(0.5) + digamma(1) + 2*log(2)
    end

    context("H() should evaluate the entropy") do
        @fact ForneyLab.H(Wishart()) --> roughly(0.7837571104739337)
    end
end

facts("Wishart converts") do
    context("Gamma should be convertible to Wishart") do
        @fact convert(Wishart, Gamma(a=1.0, b=1.0)) --> Wishart(V = [0.5].', nu = 2.0)
    end

    context("Wishart should be convertible to Gamma") do
        @fact convert(Gamma, Wishart(V = [0.5].', nu = 2.0)) --> Gamma(a=1.0, b=1.0)
        @fact_throws convert(Gamma, Wishart(V=ones(2), nu=3.0))
    end
end

