#####################
# Unit tests
#####################

facts("GaussianDistribution unit tests") do
    context("GaussianDistribution() should initialize a Gaussian distribution") do
        @fact GaussianDistribution().V => ones(1, 1)
        @fact GaussianDistribution().m => [0.0]
        @fact GaussianDistribution(m=0.0, W=1.0).W => ones(1, 1)
        @fact GaussianDistribution(xi=0.0, W=1.0).W => ones(1, 1)
        @fact typeof(GaussianDistribution(m=1.0, V=1.0).V) => Array{Float64, 2} # cast single value to matrix
        @fact_throws GaussianDistribution(V=1.0, W=1.0)
        @fact_throws GaussianDistribution(m=0.0, xi=0.0)
        @fact_throws GaussianDistribution(xi=0.0)
    end

    context("uninformative() should initialize an uninformative Gaussian distribution") do
        dist = uninformative(GaussianDistribution)
        @fact dist.m => [0.0]
        @fact dist.V => reshape([huge()],1,1)
    end

    context("Uninformative Gaussians should be equal") do
        @fact uninformative(GaussianDistribution) => uninformative(GaussianDistribution)
    end

    context("Underdetermined GaussianDistribution should be detected by isWellDefined()") do
        @fact isWellDefined(GaussianDistribution()) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, V=1.0)) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, W=1.0)) => true
        @fact isWellDefined(GaussianDistribution(xi=0.0, V=1.0)) => true
        @fact isWellDefined(GaussianDistribution(xi=0.0, W=1.0)) => true
        @fact isWellDefined(GaussianDistribution(m=0.0, xi=0.0, W=1.0, V=1.0)) => true

        dist = GaussianDistribution(m=0.0, V=1.0)
        dist.m = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, V=1.0)
        dist.V = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        dist.xi = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=0.0, W=1.0)
        dist.W = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, V=1.0, W=1.0)
        dist.m = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=0.0, xi=0.0, V=1.0)
        dist.V = nothing
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
        @fact var(GaussianDistribution(m=1.0, V=2.0)) => [2.0]
        @fact var(GaussianDistribution(m=1.0, W=2.0)) => [0.5]
        # Multivariate
        @fact mean(GaussianDistribution(m=[1.0, 2.0], V=eye(2))) => [1.0, 2.0]
        @fact mean(GaussianDistribution(xi=[1.0, 2.0], V=2.0*eye(2))) => [2.0, 4.0]
        @fact var(GaussianDistribution(m=[1.0, 2.0], V=diagm([2.0, 4.0]))) => [2.0, 4.0]
        @fact var(GaussianDistribution(m=[1.0, 2.0], W=diagm([2.0, 4.0]))) => [0.5, 0.25]
    end
    context("Inconsistent overdetermined GaussianDistribution should be detected by isConsistent()") do
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, W=1.0)) => false
        @fact isConsistent(GaussianDistribution(m=0.0, V=1.0, W=2.0)) => false
        @fact isConsistent(GaussianDistribution(m=0.0, xi=1.0, V=1.0, W=2.0)) => false
    end
end

facts("Marginal calculations for the Gaussian") do
    
    FactorGraph()

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(TerminalNode(GaussianDistribution(m=0.0, V=1.0)),
                    TerminalNode(GaussianDistribution(m=0.0, V=1.0)))
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
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
        edge = Edge(MockNode(Message(GaussianDistribution())).out, MockNode(Message(StudentsTDistribution())).out, GaussianDistribution)
        calculateMarginal!(edge)
        @fact edge.marginal => GaussianDistribution(m=0.0, W=3.0) 
    end

    context("Marginal calculation for the combination of a Gaussian and DeltaDistribution") do
        edge = Edge(MockNode(Message(GaussianDistribution())).out, MockNode(Message(DeltaDistribution(3.0))).out)
        calculateMarginal!(edge)
        @fact typeof(edge.marginal) <: DeltaDistribution => true
        @fact mean(edge.marginal) => 3.0 
    end
end

facts("GaussianDistribution converts") do
    context("DeltaDistribution should be convertible to GaussianDistribution with 0 variance") do
        @fact convert(GaussianDistribution, DeltaDistribution(3.0)) => GaussianDistribution(m=3.0, V=0.0) # Floating point number
        @fact_throws convert(GaussianDistribution, DeltaDistribution(:some_discrete_type))                # Discrete mode is not supported
        @fact_throws convert(GaussianDistribution, DeltaDistribution(ones(2,2)))                          # Matrices are not supported
        @fact convert(GaussianDistribution, DeltaDistribution([1.0, 2.0])) => GaussianDistribution(m=[1.0, 2.0], V=zeros(2,2)) # Vector
    end
end