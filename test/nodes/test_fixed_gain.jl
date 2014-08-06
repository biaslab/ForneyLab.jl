#####################
# Unit tests
#####################

facts("FixedGainNode unit tests") do
    context("FixedGainNode() should initialize a FixedGainNode with 2 interfaces") do
        node = FixedGainNode([1.0])
        @fact typeof(node) => FixedGainNode
        @fact length(node.interfaces) => 2
        @fact node.in1 => node.interfaces[1]
        @fact node.out => node.interfaces[2]
        @fact typeof(node.A) => Array{Float64, 2} # cast single value to matrix
    end

    context("FixedGainNode() should allow parameters to be clamped") do
        # Fix in1
        node = FixedGainNode(in1=GaussianDistribution())
        @fact typeof(node.in1.partner.node) => ForneyLab.ClampNode
        @fact node.in1.partner.message.payload => GaussianDistribution()
        # Fix out
        node = FixedGainNode(out=GaussianDistribution())
        @fact typeof(node.out.partner.node) => ForneyLab.ClampNode
        @fact node.out.partner.message.payload => GaussianDistribution()
    end

    context("FixedGainNode should propagate a Float") do
        # Backward message
        validateOutboundMessage(FixedGainNode([2.0]), 
                                1, 
                                Array{Float64, 2}, 
                                [nothing, Message(3.0)],
                                reshape([1.5], 1, 1))
        # Forward message
        validateOutboundMessage(FixedGainNode([2.0]), 
                                2, 
                                Array{Float64, 2}, 
                                [Message(3.0), nothing],
                                reshape([6.0], 1, 1))
    end

    context("FixedGainNode should propagate an Array") do
        # Backward message
        validateOutboundMessage(FixedGainNode([2.0]), 
                                1, 
                                Array{Float64, 1}, 
                                [nothing, Message([3.0])],
                                [1.5])
        # Forward message
        validateOutboundMessage(FixedGainNode([2.0]), 
                                2, 
                                Array{Float64, 1}, 
                                [Message([3.0]), nothing],
                                [6.0])
    end

    context("FixedGainNode should propagate a univariate GaussianDistribution") do
        # The following tests on the update rules correspond to nodes 3 and 4 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        context("Univariate GaussianDistribution with (m,V) parametrization") do
            # Backward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=3.0, V=5.0))],
                                    GaussianDistribution(m=1.5, V=1.25))
            # Forward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=3.0, V=5.0)), nothing],
                                    GaussianDistribution(m=6.0, V=20.0))
        end
        context("Univariate GaussianDistribution with (m,V=0) parametrization") do
            # Backward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=3.0, V=0.0))],
                                    GaussianDistribution(m=1.5, V=0.0))
            # Forward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=3.0, V=0.0)), nothing],
                                    GaussianDistribution(m=6.0, V=0.0))
        end
        context("Univariate GaussianDistribution with (m,W) parametrization") do
            # Backward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=3.0, W=2.0))],
                                    GaussianDistribution(m=1.5, W=8.0))
            # Forward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=3.0, W=2.0)), nothing],
                                    GaussianDistribution(m=6.0, W=0.5))
        end
        context("Univariate GaussianDistribution with (xi,W) parametrization") do
            # Backward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(xi=6.0, W=2.0))],
                                    GaussianDistribution(xi=12.0, W=8.0))
            # Forward message
            validateOutboundMessage(FixedGainNode([2.0]), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(xi=6.0, W=2.0)), nothing],
                                    GaussianDistribution(xi=3.0, W=0.5))
        end
    end

    context("FixedGainNode should propagate a multivariate GaussianDistribution") do
        # The following tests on the update rules correspond to nodes 3 and 4 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        A = reshape([   3.0, 2.0, 1.0,
                        2.0, 3.0, 2.0,
                        1.0, 2.0, 3.0], 3, 3)
        context("Multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Backward message
            validateOutboundMessage(FixedGainNode(A), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, V=variance))],
                                    GaussianDistribution(m=inv(A) * mean, V=inv(A) * variance * inv(A)'))
            # Forward message
            validateOutboundMessage(FixedGainNode(A), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, V=variance)), nothing],
                                    GaussianDistribution(m=A * mean, V=A * variance * A'))
        end
        context("Multivariate GaussianDistribution with (m,V=0) parametrization") do
            mean = [1.0:3.0]
            variance = zeros(3, 3)
            # Backward message
            validateOutboundMessage(FixedGainNode(A), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, V=variance))],
                                    GaussianDistribution(m=inv(A) * mean, V=inv(A) * variance * inv(A)'))
            # Forward message
            validateOutboundMessage(FixedGainNode(A), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, V=variance)), nothing],
                                    GaussianDistribution(m=A * mean, V=A * variance * A'))
        end
        context("Multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            # Backward message
            validateOutboundMessage(FixedGainNode(A), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, W=precision))],
                                    GaussianDistribution(m=inv(A) * mean, W=A' * precision * A))
            # Forward message
            validateOutboundMessage(FixedGainNode(A), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), nothing],
                                    GaussianDistribution(m=A * mean, W=inv(A)' * precision * inv(A)))
        end
    end
end