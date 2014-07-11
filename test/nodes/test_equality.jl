#####################
# Unit tests
#####################

facts("EqualityNode unit tests") do
   context("EqualityNode() should initialize an EqualityNode with 3 interfaces") do
        node = EqualityNode()
        @fact typeof(node) => EqualityNode
        @fact length(node.interfaces) => 3
    end

   context("EqualityNode(N) should initialize an EqualityNode with N (>2) interfaces") do
        @fact length(EqualityNode(5).interfaces) => 5
        @fact_throws EqualityNode(2)
    end
end

#####################
# Integration tests
#####################

facts("EqualityNode integration tests") do
    context("EqualityNode should propagate a univariate GaussianDistribution") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of updateNodeMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("Univariate GaussianDistribution with (m,V) parametrization") do
            inbound_dist = GaussianDistribution(m=[3.0], V=[5.0])
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            W = inv(inbound_dist.V)
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureMVParametrization!(msg.value)
            @fact isApproxEqual(msg.value.m, (pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m))) => true
            @fact isApproxEqual(msg.value.V, (inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)) => true
        end
        context("Univariate GaussianDistribution with (m,W) parametrization") do
            inbound_dist = GaussianDistribution(m=[3.0], W=[0.2])
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            W = inbound_dist.W
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureMWParametrization!(msg.value)
            @fact isApproxEqual(msg.value.m, (pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m))) => true
            @fact isApproxEqual(msg.value.W, (W + W)) => true
        end
        context("Univariate GaussianDistribution with (xi,V) parametrization") do
            inbound_dist = GaussianDistribution(xi=[0.6], V=[5.0])
            xi = inbound_dist.xi
            W = inbound_dist.W
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureXiVParametrization!(msg.value)
            @fact isApproxEqual(msg.value.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.value.V, (inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)) => true
        end
        context("Univariate GaussianDistribution with (xi,W) parametrization") do
            inbound_dist = GaussianDistribution(xi=[0.6], W=[0.2])
            xi = inbound_dist.xi
            W = inbound_dist.W
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureXiWParametrization!(msg.value)
            @fact isApproxEqual(msg.value.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.value.W, (W + W)) => true
        end
    end

    context("EqualityNode should propagate a multivariate GaussianDistribution") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of updateNodeMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("Multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = GaussianDistribution(m=mean, V=variance)
            W = inv(inbound_dist.V)
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureMVParametrization!(msg.value)
            @fact isApproxEqual(msg.value.m, (pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m))) => true
            @fact isApproxEqual(msg.value.V, (inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)) => true
        end
        context("Multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = GaussianDistribution(m=mean, W=precision)
            W = inbound_dist.W
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureMWParametrization!(msg.value)
            @fact isApproxEqual(msg.value.m, (pinv(W + W) * (W*inbound_dist.m + W*inbound_dist.m))) => true
            @fact isApproxEqual(msg.value.W, (W + W)) => true
        end
        context("Multivariate GaussianDistribution with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_dist = GaussianDistribution(xi=inv(variance)*[1.0:3.0], V=variance)
            xi = inbound_dist.xi
            V = inbound_dist.V
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureXiVParametrization!(msg.value)
            @fact isApproxEqual(msg.value.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.value.V, (inbound_dist.V * pinv(inbound_dist.V + inbound_dist.V) * inbound_dist.V)) => true
        end
        context("Multivariate GaussianDistribution with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_dist = GaussianDistribution(xi=precision*[1.0:3.0], W=precision)
            xi = inbound_dist.xi
            W = inbound_dist.W
            node = initializeEqualityNode([Message(inbound_dist), Message(inbound_dist), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            ensureXiWParametrization!(msg.value)
            @fact isApproxEqual(msg.value.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.value.W, (W + W)) => true
        end
    end

    context("EqualityNode should propagate an arbitrary message") do
        # Equality constraint node should work for arbitraty messages, although not really useful.
        # Outbound message is equal to the inbound messages if not all inbound messages are equal.
        # Otherwise, the outbound message is Message(0.0)

        # Equal scalars
        node = initializeEqualityNode([Message(1.0), Message(1.0), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Float64)
        @fact node.interfaces[3].message => msg
        @fact msg.value => 1.0
        # Unequal scalars
        node = initializeEqualityNode([Message(1.0), Message(1.1), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Float64)
        @fact msg.value => 0.0
        # Equal matrices
        node = initializeEqualityNode([Message(ones(2,2)), Message(ones(2,2)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Array{Float64})
        @fact node.interfaces[3].message => msg
        @fact msg.value => ones(2,2)
        # Unequal matrices (different values) should give zeros matrix
        node = initializeEqualityNode([Message(ones(2,2)), Message(4.0*ones(2,2)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Array{Float64})
        @fact node.interfaces[3].message => msg
        @fact msg.value => zeros(2,2)
        # Unequal matrices (different size) should give zeros matrix
        node = initializeEqualityNode([Message(ones(2,2)), Message(ones(3,3)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Array{Float64})
        @fact node.interfaces[3].message => msg
        @fact maximum(msg.value) => 0.0
    end

    context("EqualityNode should propagate a GammaMessage") do
        node = initializeEqualityNode([Message(GammaDistribution()), Message(GammaDistribution()), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GammaDistribution)
        @fact node.interfaces[3].message => msg
        @fact msg.value.a => 1.0
        @fact msg.value.b => 2.0
    end

    context("EqualityNode should propagate an InverseGammaMessage") do
        node = initializeEqualityNode([Message(InverseGammaDistribution()), Message(InverseGammaDistribution()), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, InverseGammaDistribution)
        @fact node.interfaces[3].message => msg
        @fact msg.value.a => 3.0
        @fact msg.value.b => 2.0
    end
end