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
        # In the tests, we use the exact rules from Korl. The actual implementation of calculateMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        node = EqualityNode()
        context("Univariate GaussianMessage with (m,V) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[5.0])
            W = inv(inbound_msg.V)
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureMVParametrization!(msg)
            @fact maximum(abs(msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))) < epsilon => true
            @fact maximum(abs(msg.V - (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V))) < epsilon => true
        end
        context("Univariate GaussianMessage with (m,W) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], W=[0.2])
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureMWParametrization!(msg)
            @fact maximum(abs(msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))) < epsilon => true
            @fact maximum(abs(msg.W - (W + W))) < epsilon => true
        end
        context("Univariate GaussianMessage with (xi,V) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], V=[5.0])
            xi = inbound_msg.xi
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureXiVParametrization!(msg)
            @fact maximum(abs(msg.xi - (xi + xi))) < epsilon => true
            @fact maximum(abs(msg.V - (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V))) < epsilon => true
        end
        context("Univariate GaussianMessage with (xi,W) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], W=[0.2])
            xi = inbound_msg.xi
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureXiWParametrization!(msg)
            @fact maximum(abs(msg.xi - (xi + xi))) < epsilon => true
            @fact maximum(abs(msg.W - (W + W))) < epsilon => true
        end
    end

    context("EqualityNode should propagate a multivariate GaussianMessage") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of calculateMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        node = EqualityNode()
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            W = inv(inbound_msg.V)
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureMVParametrization!(msg)
            @fact maximum(abs(msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))) < epsilon => true
            @fact maximum(abs(msg.V - (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V))) < epsilon => true
        end
        context("Multivariate GaussianMessage with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(m=mean, W=precision)
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureMWParametrization!(msg)
            @fact maximum(abs(msg.m - (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m)))) < epsilon => true
            @fact maximum(abs(msg.W - (W + W))) < epsilon => true
        end
        context("Multivariate GaussianMessage with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(xi=inv(variance)*[1.0:3.0], V=variance)
            xi = inbound_msg.xi
            V = inbound_msg.V
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureXiVParametrization!(msg)
            @fact maximum(abs(msg.xi - (xi + xi))) < epsilon => true
            @fact maximum(abs(msg.V - (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V))) < epsilon => true
        end
        context("Multivariate GaussianMessage with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(xi=precision*[1.0:3.0], W=precision)
            xi = inbound_msg.xi
            W = inbound_msg.W
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            ensureXiWParametrization!(msg)
            @fact maximum(abs(msg.xi - (xi + xi))) < epsilon => true
            @fact maximum(abs(msg.W - (W + W))) < epsilon => true
        end
    end
end