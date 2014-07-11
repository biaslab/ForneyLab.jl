#####################
# Unit tests
#####################

facts("GaussianDistribution unit tests") do
    context("GaussianDistribution() should initialize a Gaussian distribution") do
        @fact GaussianDistribution().V => ones(1, 1)
        @fact GaussianDistribution().m => [0.0]
        @fact GaussianDistribution(m=[0.0], W=[1.0]).W => ones(1, 1)
        @fact GaussianDistribution(xi=[0.0], W=[1.0]).W => ones(1, 1)
        @fact typeof(GaussianDistribution(m=[1.0], V=[1.0]).V) => Array{Float64, 2} # cast single value to matrix
        @fact_throws GaussianDistribution(V=[1.0], W=[1.0])
        @fact_throws GaussianDistribution(m=[0.0], xi=[0.0])
        @fact_throws GaussianDistribution(xi=[0.0])
    end

    context("uninformative() should initialize an uninformative Gaussian distribution") do
        dist = uninformative(GaussianDistribution)
        @fact dist.m => [0.0]
        @fact dist.V => reshape([1000.0],1,1)
    end

    context("Underdetermined GaussianDistribution should be detected by isWellDefined()") do
        @fact isWellDefined(GaussianDistribution()) => true
        @fact isWellDefined(GaussianDistribution(m=[0.0], V=[1.0])) => true
        @fact isWellDefined(GaussianDistribution(m=[0.0], W=[1.0])) => true
        @fact isWellDefined(GaussianDistribution(xi=[0.0], V=[1.0])) => true
        @fact isWellDefined(GaussianDistribution(xi=[0.0], W=[1.0])) => true
        @fact isWellDefined(GaussianDistribution(m=[0.0], xi=[0.0], W=[1.0], V=[1.0])) => true

        dist = GaussianDistribution(m=[0.0], V=[1.0])
        dist.m = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=[0.0], V=[1.0])
        dist.V = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=[0.0], W=[1.0])
        dist.xi = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(xi=[0.0], W=[1.0])
        dist.W = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=[0.0], V=[1.0], W=[1.0])
        dist.m = nothing
        @fact isWellDefined(dist) => false

        dist = GaussianDistribution(m=[0.0], xi=[0.0], V=[1.0])
        dist.V = nothing
        @fact isWellDefined(dist) => false
    end
    context("Conversions between valid parametrizations of a GaussianDistribution should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(m=[0.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(m=[0.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(m=[0.0], V=[1.0]))) => true
        # Defined as (m,W)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(m=[0.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(m=[0.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(m=[0.0], W=[1.0]))) => true
        # Defined as (xi,V)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(xi=[2.0], V=[1.0]))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(xi=[2.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianDistribution(xi=[2.0], V=[1.0]))) => true
        # Defined as (xi,W)
        @fact isConsistent(ensureMVParametrization!(GaussianDistribution(xi=[2.0], W=[1.0]))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianDistribution(xi=[2.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianDistribution(xi=[2.0], W=[1.0]))) => true
    end
    context("mean(GaussianDistribution) and var(GaussianDistribution) should return correct result") do
        # Univariate
        @fact mean(GaussianDistribution(m=[1.0], V=[1.0])) => [1.0]
        @fact mean(GaussianDistribution(xi=[1.0], V=[2.0])) => [2.0]
        @fact var(GaussianDistribution(m=[1.0], V=[2.0])) => [2.0]
        @fact var(GaussianDistribution(m=[1.0], W=[2.0])) => [0.5]
        # Multivariate
        @fact mean(GaussianDistribution(m=[1.0, 2.0], V=eye(2))) => [1.0, 2.0]
        @fact mean(GaussianDistribution(xi=[1.0, 2.0], V=2.0*eye(2))) => [2.0, 4.0]
        @fact var(GaussianDistribution(m=[1.0, 2.0], V=diagm([2.0, 4.0]))) => [2.0, 4.0]
        @fact var(GaussianDistribution(m=[1.0, 2.0], W=diagm([2.0, 4.0]))) => [0.5, 0.25]
    end
    context("Inconsistent overdetermined GaussianDistribution should be detected by isConsistent()") do
        @fact isConsistent(GaussianDistribution(m=[0.0], xi=[1.0], W=[1.0])) => false
        @fact isConsistent(GaussianDistribution(m=[0.0], V=[1.0], W=[2.0])) => false
        @fact isConsistent(GaussianDistribution(m=[0.0], xi=[1.0], V=[1.0], W=[2.0])) => false
    end
end

facts("Marginal calculations for the Gaussian") do
    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(ConstantNode(GaussianDistribution(m=[0.0], V=[1.0])),
                    ConstantNode(GaussianDistribution(m=[0.0], V=[1.0])))
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
                                GaussianDistribution(m=[0.0], V=[1.0]),
                                GaussianDistribution(m=[0.0], V=[1.0]))
        ensureMVParametrization!(marginal_dist)
        @fact marginal_dist.m => [0.0]
        @fact isApproxEqual(marginal_dist.V, reshape([0.5], 1, 1)) => true
    end
end