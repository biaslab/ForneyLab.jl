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
        # In the tests, we use the exact rules from Korl. The actual implementation of updateNodeMessage!() will calculate
        # the (xi,W) parametrizations of the inbound messages, such that only the W and xi update rules are used in practice.
        node = EqualityNode()
        context("Univariate GaussianMessage with (m,V) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[5.0])
            W = inv(inbound_msg.V)
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            @fact node.interfaces[3].message => msg
            ensureMVParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Univariate GaussianMessage with (m,W) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], W=[0.2])
            W = inbound_msg.W
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            @fact node.interfaces[3].message => msg
            ensureMWParametrization!(msg)
            @fact isApproxEqual(msg.m, (pinv(W + W) * (W*inbound_msg.m + W*inbound_msg.m))) => true
            @fact isApproxEqual(msg.W, (W + W)) => true
        end
        context("Univariate GaussianMessage with (xi,V) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], V=[5.0])
            xi = inbound_msg.xi
            W = inbound_msg.W
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
            @fact node.interfaces[3].message => msg
            ensureXiVParametrization!(msg)
            @fact isApproxEqual(msg.xi, (xi + xi)) => true
            @fact isApproxEqual(msg.V, (inbound_msg.V * pinv(inbound_msg.V + inbound_msg.V) * inbound_msg.V)) => true
        end
        context("Univariate GaussianMessage with (xi,W) parametrization") do
            inbound_msg = GaussianMessage(xi=[0.6], W=[0.2])
            xi = inbound_msg.xi
            W = inbound_msg.W
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
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
        node = EqualityNode()
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            W = inv(inbound_msg.V)
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
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
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
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
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
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
            inbound_messages = Array(GaussianMessage, 3)
            inbound_messages[1] = inbound_msg
            inbound_messages[2] = inbound_msg
            msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
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
        node = EqualityNode()
        inbound_messages = Array(GeneralMessage, 3)
        # Equal scalars
        inbound_messages[1] = GeneralMessage(1.0)
        inbound_messages[2] = GeneralMessage(1.0)
        msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
        @fact node.interfaces[3].message => msg
        @fact msg.value => 1.0
        # Unequal scalars
        inbound_messages[2] = GeneralMessage(1.1)
        msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
        @fact msg.value => 0.0
        # Equal matrices
        inbound_messages[1] = GeneralMessage(ones(2,2))
        inbound_messages[2] = GeneralMessage(ones(2,2))
        msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
        @fact node.interfaces[3].message => msg
        @fact msg.value => ones(2,2)
        # Unequal matrices (different values) should give zeros matrix
        inbound_messages[2] = GeneralMessage(4.0*ones(2,2))
        msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
        @fact node.interfaces[3].message => msg
        @fact msg.value => zeros(2,2)
        # Unequal matrices (different size) should give zeros matrix
        inbound_messages[2] = GeneralMessage(ones(3,3))
        msg = ForneyLab.updateNodeMessage!(3, node, inbound_messages)
        @fact node.interfaces[3].message => msg
        @fact maximum(msg.value) => 0.0
    end
end