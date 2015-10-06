#####################
# Unit tests
#####################

facts("GaussianDistribution unit tests") do
    context("GaussianDistribution() should initialize a Gaussian distribution") do
        @fact GaussianDistribution().V => eye(1)
        @fact GaussianDistribution().m => [0.0]
        @fact GaussianDistribution(m=0.0, W=1.0).W => ones(1, 1)
        @fact GaussianDistribution(xi=0.0, W=1.0).W => ones(1, 1)
        @fact typeof(GaussianDistribution(m=1.0, V=1.0).V) => Array{Float64, 2} # cast single value to matrix
        @fact GaussianDistribution(m=[0.0, 0.0], V=eye(2)).V => eye(2) # multivariate
        @fact_throws GaussianDistribution(V=1.0, W=1.0)
        @fact_throws GaussianDistribution(m=0.0, xi=0.0)
        @fact_throws GaussianDistribution(xi=0.0)
        @fact_throws GaussianDistribution(m=[0.0, 0.0], V=[0. 0.;0. 0.])  # V should be positive definite
        @fact_throws GaussianDistribution(m=[0.0, 0.0], V=[Inf 0.;0. 1.]) # V should be bounded
        @fact_throws GaussianDistribution(m=[0.0, 0.0], W=[0. 0.;0. 0.])  # W should be positive definite
        @fact_throws GaussianDistribution(m=[0.0, 0.0], W=[Inf 0.;0. 1.]) # W should be bounded
    end

    context("vague() should initialize a vague (almost uninformative) Gaussian distribution") do
        dist = vague(GaussianDistribution)
        @fact dist.m => [0.0]
        @fact dist.V => reshape([huge],1,1)
    end

    context("isProper() should indicate whether a Gaussian distribution is proper") do
        @fact isProper(GaussianDistribution()) => true
        @fact isProper(GaussianDistribution(m=0.0, V=-1.0)) => false
        @fact isProper(GaussianDistribution(m=0.0, W=-1.0)) => false
    end

    context("Vague Gaussians should be equal") do
        @fact vague(GaussianDistribution) => vague(GaussianDistribution)
    end

    context("GaussianDistribution can be sampled") do
        @fact typeof(sample(GaussianDistribution(m=[1.2, 2.7], V=[2.0 -0.5; -0.5 1.5]))) => Array{Float64, 1}
    end

    context("Underdetermined GaussianDistribution should be detected by isWellDefined()") do
        @fact isWellDefined(GaussianDistribution()) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, V=1.0)) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, W=1.0)) => true
        @fact isWellDefined(GaussianDistribution(xi=0.0, V=1.0)) => true
        @fact isWellDefined(GaussianDistribution(xi=0.0, W=1.0)) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, xi=0.0, W=1.0, V=1.0)) => true

        dist = GaussianDistribution(m=0.0, V=1.0)
        invalidate!(dist.m)
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, V=1.0)
        invalidate!(dist.V)
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        invalidate!(dist.xi)
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        invalidate!(dist.W)
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, V=1.0, W=1.0)
        invalidate!(dist.m)
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, xi=0.0, V=1.0)
        invalidate!(dist.V)
        @fact isWellDefined(dist) => false
    end
    context("Conversions between valid parametrizations of a GaussianDistribution should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(m=0.0, V=1.0))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(m=0.0, V=1.0))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(m=0.0, V=1.0))) => true
        # Defined as (m,W)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(m=0.0, W=1.0))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(m=0.0, W=1.0))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(m=0.0, W=1.0))) => true
        # Defined as (xi,V)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(xi=2.0, V=1.0))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(xi=2.0, V=1.0))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(xi=2.0, V=1.0))) => true
        # Defined as (xi,W)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(xi=2.0, W=1.0))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(xi=2.0, W=1.0))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(xi=2.0, W=1.0))) => true
    end
    context("mean(GaussianDistribution) and var(GaussianDistribution) should return correct result") do
        # Univariate
        @fact mean(GaussianDistribution(m=1.0, V=1.0)) => [1.0]
        @fact mean(GaussianDistribution(xi=1.0, V=2.0)) => [2.0]
        @fact isValid(mean(GaussianDistribution(xi=1.0, V=-2.0))) => false
        @fact var(GaussianDistribution(m=1.0, V=2.0)) => [2.0]
        @fact var(GaussianDistribution(m=1.0, W=2.0)) => [0.5]
        @fact isValid(var(GaussianDistribution(m=1.0, W=-2.0))) => false
        # Multivariate
        @fact mean(GaussianDistribution(m=[1.0, 2.0], V=eye(2))) => [1.0, 2.0]
        @fact mean(GaussianDistribution(xi=[1.0, 2.0], V=2.0*eye(2))) => [2.0, 4.0]
        @fact isValid(mean(GaussianDistribution(xi=[1.0, 2.0], V=-2.0*eye(2)))) => false
        @fact var(GaussianDistribution(m=[1.0, 2.0], V=diagm([2.0, 4.0]))) => [2.0, 4.0]
        @fact var(GaussianDistribution(m=[1.0, 2.0], W=diagm([2.0, 4.0]))) => [0.5, 0.25]
        @fact isValid(var(GaussianDistribution(m=[1.0, 2.0], W=diagm([-2.0, 4.0])))) => false
    end
    context("Inconsistent overdetermined GaussianDistribution should be detected by isConsistent()") do
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, W=1.0)) => false
        @fact isConsistent(GaussianDistribution(m=0.0, V=1.0, W=2.0)) => false
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, V=1.0, W=2.0)) => false
    end
end

facts("Marginal calculations for the Gaussian") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        initializePairOfTerminalNodes(GaussianDistribution(m=0.0, V=1.0), GaussianDistribution(m=0.0, V=1.0))
        edge = n(:t1).i[:out].edge
        n(:t1).i[:out].message = Message(GaussianDistribution(m=0.0, V=1.0))
        n(:t2).i[:out].message = Message(GaussianDistribution(m=0.0, V=1.0))
        marginal_dist = calculateMarginal!(edge)
        @fact edge.marginal => marginal_dist
        ensureMVParametrization!(marginal_dist)
        @fact edge.marginal.m => [0.0]
        @fact isApproxEqual(edge.marginal.V, reshape([0.5], 1, 1)) => true
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_dist = calculateMarginal(
                                GaussianDistribution(m=0.0, V=1.0),
                                GaussianDistribution(m=0.0, V=1.0))
        ensureMVParametrization!(marginal_dist)
        @fact marginal_dist.m => [0.0]
        @fact isApproxEqual(marginal_dist.V, reshape([0.5], 1, 1)) => true
    end

    context("Marginal calculation for the combination of a Gaussian and student's t-distribution") do
        initializePairOfTerminalNodes(GaussianDistribution(), StudentsTDistribution(m=1.0, lambda=2.0, nu=4.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal => GaussianDistribution(m=0.5, W=2.0) 
    end

    context("Marginal calculation for the combination of a Gaussian and DeltaDistribution") do
        initializePairOfTerminalNodes(GaussianDistribution(), DeltaDistribution(3.0))
        edge = n(:t1).i[:out].edge
        calculateMarginal!(edge)
        @fact edge.marginal => DeltaDistribution(3.0)
    end
end

facts("GaussianDistribution converts") do
    context("DeltaDistribution should be convertible to GaussianDistribution with negligible variance") do
        @fact convert(GaussianDistribution, DeltaDistribution(3.0)) => GaussianDistribution(m=3.0, V=tiny) # Floating point number
        @fact convert(GaussianDistribution, DeltaDistribution([1.0, 2.0])) => GaussianDistribution(m=[1.0, 2.0], V=tiny*eye(2)) # Vector
    end

    context("Message{DeltaDistribution} should be convertible to Message{GaussianDistribution} with negligible variance") do
        @fact convert(Message{GaussianDistribution}, Message(DeltaDistribution(3.0))) => Message(GaussianDistribution(m=3.0, V=tiny)) # Floating point number
        @fact convert(Message{GaussianDistribution}, Message(DeltaDistribution([1.0, 2.0]))) => Message(GaussianDistribution(m=[1.0, 2.0], V=tiny*eye(2))) # Vector
    end
end
