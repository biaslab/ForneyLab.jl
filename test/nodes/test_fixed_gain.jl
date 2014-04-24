facts("FixedGainNode") do
    context("FixedGainNode() should initialize a FixedGainNode with 2 interfaces") do
        node = FixedGainNode([1.0])
        @fact typeof(node) => FixedGainNode
        @fact length(node.interfaces) => 2
        @fact node.in1 => node.interfaces[1]
        @fact node.out => node.interfaces[2]
        @fact typeof(node.A) => Array{Float64, 2} # cast single value to matrix
    end

    context("FixedGainNode should propagate a GeneralMessage") do
        node = FixedGainNode([2.0])
        inbound_msg = GeneralMessage([3.0])
        # Backward message
        msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => [1.5]
        # Forward message
        msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
        @fact node.interfaces[2].message => msg
        @fact node.interfaces[2].message.value => [6]
    end

    context("FixedGainNode should propagate a univariate GaussianMessage") do
        node = FixedGainNode([2.0])
        context("Univariate GaussianMessage with (m,V) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[5.0])
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [1.5]
            @fact node.interfaces[1].message.V => reshape([1.25], 1, 1)
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [6.0]
            @fact node.interfaces[2].message.V => reshape([20.0], 1, 1)
        end
        context("Univariate GaussianMessage with (m,V=0) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], V=[0.0])
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [1.5]
            @fact node.interfaces[1].message.V => reshape([0.0], 1, 1)
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [6.0]
            @fact node.interfaces[2].message.V => reshape([0.0], 1, 1)
        end
        context("Univariate GaussianMessage with (m,W) parametrization") do
            inbound_msg = GaussianMessage(m=[3.0], W=[2.0])
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [1.5]
            @fact node.interfaces[1].message.W => reshape([8.0], 1, 1)
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [6.0]
            @fact node.interfaces[2].message.W => reshape([0.5], 1, 1)
        end
        context("Univariate GaussianMessage with (xi,W) parametrization") do
            inbound_msg = GaussianMessage(xi=[6.0], W=[2.0])
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.xi => [12.0]
            @fact node.interfaces[1].message.W => reshape([8.0], 1, 1)
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.xi => [3.0]
            @fact node.interfaces[2].message.W => reshape([0.5], 1, 1)
        end
    end

    context("FixedGainNode should propagate a multivariate GaussianMessage") do
        A = reshape([   3.0, 2.0, 1.0,
                        2.0, 3.0, 2.0,
                        1.0, 2.0, 3.0], 3, 3)
        node = FixedGainNode(A)
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => inv(A) * mean
            @fact node.interfaces[1].message.V => inv(A) * variance * inv(A)'
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => A * mean
            @fact node.interfaces[2].message.V => A * variance * A'
        end
        context("Multivariate GaussianMessage with (m,V=0) parametrization") do
            mean = [1.0:3.0]
            variance = zeros(3, 3)
            inbound_msg = GaussianMessage(m=mean, V=variance)
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => inv(A) * mean
            @fact node.interfaces[1].message.V => inv(A) * variance * inv(A)'
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => A * mean
            @fact node.interfaces[2].message.V => A * variance * A'
        end
        context("Multivariate GaussianMessage with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            inbound_msg = GaussianMessage(m=mean, W=precision)
            # Backward message
            msg = calculateMessage!(1, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => inv(A) * mean
            @fact node.interfaces[1].message.W => A' * precision * A
            # Forward message
            msg = calculateMessage!(2, node, [inbound_msg, inbound_msg])
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => A * mean
            @fact node.interfaces[2].message.W => inv(A)' * precision * inv(A)
        end
    end
end