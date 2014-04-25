facts("EqualityNode") do
   context("EqualityNode() should initialize an EqualityNode with 3 interfaces") do
        node = EqualityNode()
        @fact typeof(node) => EqualityNode
        @fact length(node.interfaces) => 3
    end

   context("EqualityNode(N) should initialize an EqualityNode with N (>2) interfaces") do
        @fact length(EqualityNode(5).interfaces) => 5
        @fact_throws EqualityNode(2)
    end

    context("EqualityNode should propagate a univariate GaussianMessage") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        node = EqualityNode()
        delta = 2.0*eps(Float64)
        context("Univariate GaussianMessage with (m,V) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[5.0])
            W = inv(inbound_msg.V)
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact (msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))[1] < delta => true
            @fact (msg.V - (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V))[1] < delta => true
        end
        context("Univariate GaussianMessage with (m,W) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], W=[0.2])
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact (msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))[1] < delta => true
            @fact (msg.W - (W + W))[1] < delta => true
        end
        context("Univariate GaussianMessage with (xi,W) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], W=[0.2])
            xi = inbound_msg.xi
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact (msg.xi - (xi + xi))[1] < delta => true
            @fact (msg.W - (W + W))[1] < delta => true
        end
    end

    context("EqualityNode should propagate a multivariate GaussianMessage") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        node = EqualityNode()
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            # TODO
        end
        context("Multivariate GaussianMessage with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(m=mean, W=precision)
            # TODO
        end
        context("Multivariate GaussianMessage with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            weighted_mean = precision * [1.0:3.0]
            inbound_msg = GaussianMessage(xi=weighted_mean, W=precision)
            # TODO
        end
    end
end