#####################
# Unit tests
#####################

facts("GaussianDistribution unit tests") do
    context("GaussianDistribution() should initialize a Gaussian distribution") do
        @fact GaussianDistribution().V --> 1.0
        @fact GaussianDistribution().m --> 0.0
        @fact GaussianDistribution(m=0.0, W=1.0).W --> 1.0
        @fact GaussianDistribution(xi=0.0, W=1.0).W --> 1.0
        @fact_throws GaussianDistribution(V=1.0, W=1.0)
        @fact_throws GaussianDistribution(m=0.0, xi=0.0)
        @fact_throws GaussianDistribution(xi=0.0)
        @fact_throws GaussianDistribution(m=0.0, V=0.0)  # V should be positive definite
        @fact_throws GaussianDistribution(m=0.0, V=Inf) # V should be bounded
        @fact_throws GaussianDistribution(m=0.0, W=0.0)  # W should be positive definite
        @fact_throws GaussianDistribution(m=0.0, W=Inf) # W should be bounded
    end

    context("vague() should initialize a vague (almost uninformative) Gaussian distribution") do
        dist = vague(GaussianDistribution)
        @fact dist.m --> 0.0
        @fact dist.V --> huge
    end

    context("isProper() should indicate whether a Gaussian distribution is proper") do
        @fact isProper(GaussianDistribution()) --> true
        @fact isProper(GaussianDistribution(m=0.0, V=-1.0)) --> false
        @fact isProper(GaussianDistribution(m=0.0, W=-1.0)) --> false
    end

    context("GaussianDistribution can be sampled") do
        @fact typeof(sample(GaussianDistribution(m=1.0, V=2.0))) --> Float64
    end

    context("prod!() should produce correct output") do
        p = GaussianDistribution(m=0.0, V=1.0) * GaussianDistribution(m=0.0, V=1.0)
        ensureParameters!(p, (:m, :V))
        @fact p.m --> roughly(0.0)
        @fact p.V --> roughly(0.5)
        p = GaussianDistribution(xi=0.0, W=1.0) * GaussianDistribution(xi=0.0, W=1.0)
        ensureParameters!(p, (:m, :V))
        @fact p.m --> roughly(0.0)
        @fact p.V --> roughly(0.5)
        @fact GaussianDistribution() * DeltaDistribution(2.0) --> DeltaDistribution(2.0)
        @fact DeltaDistribution(2.0) *  GaussianDistribution() --> DeltaDistribution(2.0)
        @fact ForneyLab.prod!(GaussianDistribution(), DeltaDistribution(2.0), GaussianDistribution()) --> GaussianDistribution(m=2.0, V=tiny)
        @fact ForneyLab.prod!(DeltaDistribution(2.0), GaussianDistribution(), GaussianDistribution()) --> GaussianDistribution(m=2.0, V=tiny)
    end

    context("Underdetermined GaussianDistribution should be detected by isWellDefined()") do
        @fact isWellDefined(GaussianDistribution()) --> true
        @fact isWellDefined(GaussianDistribution(m=0.0, V=1.0)) --> true
        @fact isWellDefined(GaussianDistribution(m=0.0, W=1.0)) --> true
        @fact isWellDefined(GaussianDistribution(xi=0.0, V=1.0)) --> true
        @fact isWellDefined(GaussianDistribution(xi=0.0, W=1.0)) --> true
        @fact isWellDefined(GaussianDistribution(m=0.0, xi=0.0, W=1.0, V=1.0)) --> true

        dist = GaussianDistribution(m=0.0, V=1.0)
        dist.m = NaN # Invalidate
        @fact isWellDefined(dist) --> false

        dist = GaussianDistribution(m=0.0, V=1.0)
        dist.V = NaN
        @fact isWellDefined(dist) --> false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        dist.xi = NaN
        @fact isWellDefined(dist) --> false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        dist.W = NaN
        @fact isWellDefined(dist) --> false

        dist = GaussianDistribution(m=0.0, V=1.0, W=1.0)
        dist.m = NaN
        @fact isWellDefined(dist) --> false

        dist = GaussianDistribution(m=0.0, xi=0.0, V=1.0)
        dist.V = NaN
        @fact isWellDefined(dist) --> false
    end
    context("Conversions between valid parametrizations of a GaussianDistribution should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, V=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, V=1.0), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, V=1.0), (:xi, :W))) --> true
        # Defined as (m,W)
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, W=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, W=1.0), (:xi, :V))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(m=0.0, W=1.0), (:xi, :W))) --> true
        # Defined as (xi,V)
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, V=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, V=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, V=1.0), (:xi, :W))) --> true
        # Defined as (xi,W)
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, W=1.0), (:m, :V))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, W=1.0), (:m, :W))) --> true
        @fact isConsistent(ensureParameters!(GaussianDistribution(xi=2.0, W=1.0), (:xi, :V))) --> true
    end
    context("mean(GaussianDistribution) and var(GaussianDistribution) should return correct result") do
        @fact mean(GaussianDistribution(m=1.0, V=1.0)) --> 1.0
        @fact mean(GaussianDistribution(xi=1.0, V=2.0)) --> 2.0
        @fact isnan(mean(GaussianDistribution(xi=1.0, V=-2.0))) --> true
        @fact var(GaussianDistribution(m=1.0, V=2.0)) --> 2.0
        @fact var(GaussianDistribution(m=1.0, W=2.0)) --> 0.5
        @fact isnan(var(GaussianDistribution(m=1.0, W=-2.0))) --> true
    end
    context("Inconsistent overdetermined GaussianDistribution should be detected by isConsistent()") do
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, W=1.0)) --> false
        @fact isConsistent(GaussianDistribution(m=0.0, V=1.0, W=2.0)) --> false
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, V=1.0, W=2.0)) --> false
    end
end

facts("GaussianDistribution converts") do
    context("DeltaDistribution should be convertible to GaussianDistribution with negligible variance") do
        @fact convert(GaussianDistribution, DeltaDistribution(3.0)) --> GaussianDistribution(m=3.0, V=tiny) # Floating point number
    end

    context("Message{DeltaDistribution} should be convertible to Message{GaussianDistribution} with negligible variance") do
        @fact convert(Message{GaussianDistribution}, Message(DeltaDistribution(3.0))) --> Message(GaussianDistribution(m=3.0, V=tiny)) # Floating point number
    end
end
