#####################
# Unit tests
#####################

facts("Gaussian unit tests") do
    context("Gaussian() should initialize a Gaussian distribution") do
        @fact Gaussian().V --> 1.0
        @fact Gaussian().m --> 0.0
        @fact Gaussian(m=0.0, W=1.0).W --> 1.0
        @fact Gaussian(xi=0.0, W=1.0).W --> 1.0
        @fact_throws Gaussian(V=1.0, W=1.0)
        @fact_throws Gaussian(m=0.0, xi=0.0)
        @fact_throws Gaussian(xi=0.0)
        @fact_throws Gaussian(m=0.0, V=0.0)  # V should be positive definite
        @fact_throws Gaussian(m=0.0, V=Inf) # V should be bounded
        @fact_throws Gaussian(m=0.0, W=0.0)  # W should be positive definite
        @fact_throws Gaussian(m=0.0, W=Inf) # W should be bounded
        @fact pdf(Gaussian(m=1.0, V=4.0), 1.5) --> roughly(0.19333405840142462, atol=1e-8)
    end

    context("vague() should initialize a vague (almost uninformative) Gaussian distribution") do
        dist = vague(Gaussian)
        @fact dist.m --> 0.0
        @fact dist.V --> huge
    end

    context("isProper() should indicate whether a Gaussian distribution is proper") do
        @fact isProper(Gaussian()) --> true
        @fact isProper(Gaussian(m=0.0, V=-1.0)) --> false
        @fact isProper(Gaussian(m=0.0, W=-1.0)) --> false
    end

    context("Gaussian can be sampled") do
        @fact typeof(sample(Gaussian(m=1.0, V=2.0))) --> Float64
    end

    context("prod!() should produce correct output") do
        p = Gaussian(m=0.0, V=1.0) * Gaussian(m=0.0, V=1.0)
        ensureParameters!(p, (:m, :V))
        @fact p.m --> roughly(0.0)
        @fact p.V --> roughly(0.5)
        p = Gaussian(xi=0.0, W=1.0) * Gaussian(xi=0.0, W=1.0)
        ensureParameters!(p, (:m, :V))
        @fact p.m --> roughly(0.0)
        @fact p.V --> roughly(0.5)
        @fact Gaussian() * Delta(2.0) --> Delta(2.0)
        @fact Delta(2.0) *  Gaussian() --> Delta(2.0)
        @fact ForneyLab.prod!(Gaussian(), Delta(2.0), Gaussian()) --> Gaussian(m=2.0, V=tiny)
        @fact ForneyLab.prod!(Delta(2.0), Gaussian(), Gaussian()) --> Gaussian(m=2.0, V=tiny)
    end

    context("Underdetermined Gaussian should be detected by isWellDefined()") do
        @fact isWellDefined(Gaussian()) --> true
        @fact isWellDefined(Gaussian(m=0.0, V=1.0)) --> true
        @fact isWellDefined(Gaussian(m=0.0, W=1.0)) --> true
        @fact isWellDefined(Gaussian(xi=0.0, V=1.0)) --> true
        @fact isWellDefined(Gaussian(xi=0.0, W=1.0)) --> true
        @fact isWellDefined(Gaussian(m=0.0, xi=0.0, W=1.0, V=1.0)) --> true

        dist = Gaussian(m=0.0, V=1.0)
        dist.m = NaN # Invalidate
        @fact isWellDefined(dist) --> false

        dist = Gaussian(m=0.0, V=1.0)
        dist.V = NaN
        @fact isWellDefined(dist) --> false

        dist = Gaussian(xi=0.0, W=1.0)
        dist.xi = NaN
        @fact isWellDefined(dist) --> false

        dist = Gaussian(xi=0.0, W=1.0)
        dist.W = NaN
        @fact isWellDefined(dist) --> false

        dist = Gaussian(m=0.0, V=1.0, W=1.0)
        dist.m = NaN
        @fact isWellDefined(dist) --> false

        dist = Gaussian(m=0.0, xi=0.0, V=1.0)
        dist.V = NaN
        @fact isWellDefined(dist) --> false
    end
    context("Conversions between valid parametrizations of a Gaussian should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, V=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, V=1.0), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, V=1.0), (:xi, :W))) --> true
        # Defined as (m,W)
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, W=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, W=1.0), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(m=0.0, W=1.0), (:xi, :W))) --> true
        # Defined as (xi,V)
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, V=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, V=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, V=1.0), (:xi, :W))) --> true
        # Defined as (xi,W)
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, W=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, W=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(Gaussian(xi=2.0, W=1.0), (:xi, :V))) --> true
    end
    context("mean(Gaussian) and var(Gaussian) should return correct result") do
        @fact mean(Gaussian(m=1.0, V=1.0)) --> 1.0
        @fact mean(Gaussian(xi=1.0, V=2.0)) --> 2.0
        @fact isnan(mean(Gaussian(xi=1.0, V=-2.0))) --> true
        @fact var(Gaussian(m=1.0, V=2.0)) --> 2.0
        @fact var(Gaussian(m=1.0, W=2.0)) --> 0.5
        @fact isnan(var(Gaussian(m=1.0, W=-2.0))) --> true
    end
    context("Inconsistent overdetermined Gaussian should be detected by isConsistent()") do
        @fact isConsistent(Gaussian(m=0.0, xi=1.0, W=1.0)) --> false
        @fact isConsistent(Gaussian(m=0.0, V=1.0, W=2.0)) --> false
        @fact isConsistent(Gaussian(m=0.0, xi=1.0, V=1.0, W=2.0)) --> false
    end
end

facts("Gaussian converts") do
    context("Delta should be convertible to Gaussian with negligible variance") do
        @fact convert(Gaussian, Delta(3.0)) --> Gaussian(m=3.0, V=tiny) # Floating point number
    end

    context("Message{Delta} should be convertible to Message{Gaussian} with negligible variance") do
        @fact convert(Message{Gaussian}, Message(Delta(3.0))) --> Message(Gaussian(m=3.0, V=tiny)) # Floating point number
    end
end
