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
    context("EqualityNode should propagate a univariate GaussianMessage") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of updateNodeMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("Univariate GaussianMessage with (m,V) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[5.0])
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            W = inv(inbound_msg.V)
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureMVParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Univariate GaussianMessage with (m,W) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], W=[0.2])
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            W = inbound_msg.W
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureMWParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.W, (W + W)) => true
        end
        context("Univariate GaussianMessage with (xi,V) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], V=[5.0])
            xi = inbound_msg.xi
            W = inbound_msg.W
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureXiVParametrization!(msg)
            @fact isApproxEqual(msg.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Univariate GaussianMessage with (xi,W) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], W=[0.2])
            xi = inbound_msg.xi
            W = inbound_msg.W
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureXiWParametrization!(msg)
            @fact isApproxEqual(msg.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.W, (W + W)) => true
        end
    end

    context("EqualityNode should propagate a multivariate GaussianMessage") do
        # The following tests on the update rules correspond to node 1 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        # In the tests, we use the exact rules from Korl. The actual implementation of updateNodeMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            W = inv(inbound_msg.V)
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureMVParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Multivariate GaussianMessage with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(m=mean, W=precision)
            W = inbound_msg.W
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureMWParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.W, (W + W)) => true
        end
        context("Multivariate GaussianMessage with (xi,V) parametrization") do
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(xi=inv(variance)*[1.0:3.0], V=variance)
            xi = inbound_msg.xi
            V = inbound_msg.V
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureXiVParametrization!(msg)
            @fact isApproxEqual(msg.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Multivariate GaussianMessage with (xi,W) parametrization") do
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(xi=precision*[1.0:3.0], W=precision)
            xi = inbound_msg.xi
            W = inbound_msg.W
            node = initializeEqualityNode([inbound_msg, inbound_msg, nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            ensureXiWParametrization!(msg)
            @fact isApproxEqual(msg.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.W, (W + W)) => true
        end
    end

    context("EqualityNode should propagate a GeneralMessage") do
        # Equality constraint node should work for GeneralMessages, although not really useful.
        # Outbound message is equal to the inbound messages if not all inbound messages are equal.
        # Otherwise, the outbound message is GeneralMessage(0.0)

        # Equal scalars
        node = initializeEqualityNode([GeneralMessage(1.0), GeneralMessage(1.0), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact node.interfaces[3].message => msg
        @fact msg.value => 1.0
        # Unequal scalars
        node = initializeEqualityNode([GeneralMessage(1.0), GeneralMessage(1.1), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact msg.value => 0.0
        # Equal matrices
        node = initializeEqualityNode([GeneralMessage(ones(2,2)), GeneralMessage(ones(2,2)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact node.interfaces[3].message => msg
        @fact msg.value => ones(2,2)
        # Unequal matrices (different values) should give zeros matrix
        node = initializeEqualityNode([GeneralMessage(ones(2,2)), GeneralMessage(4.0*ones(2,2)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact node.interfaces[3].message => msg
        @fact msg.value => zeros(2,2)
        # Unequal matrices (different size) should give zeros matrix
        node = initializeEqualityNode([GeneralMessage(ones(2,2)), GeneralMessage(ones(3,3)), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact node.interfaces[3].message => msg
        @fact maximum(msg.value) => 0.0
    end

    context("EqualityNode should propagate a GammaMessage") do
        node = initializeEqualityNode([GammaMessage(inverted=true), GammaMessage(inverted=true), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GammaMessage)
        @fact node.interfaces[3].message => msg
        @fact msg.a => 3.0
        @fact msg.b => 2.0
        @fact msg.inverted => true
    end
end