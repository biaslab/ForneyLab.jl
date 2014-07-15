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
end

#####################
# Integration tests
#####################

facts("FixedGainNode integration tests") do
    context("FixedGainNode should propagate a Float") do
        # Backward message
        node = initializeFixedGainNode([2.0], [nothing, Message(3.0)])
        msg = ForneyLab.updateNodeMessage!(1, node, Array{Float64})
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => reshape([1.5], 1, 1)
        # Forward message
        node = initializeFixedGainNode([2.0], [Message(3.0), nothing])
        msg = ForneyLab.updateNodeMessage!(2, node, Array{Float64})
        @fact node.interfaces[2].message => msg
        @fact node.interfaces[2].message.value => reshape([6.0], 1, 1)
    end

    context("FixedGainNode should propagate an Array") do
        # Backward message
        node = initializeFixedGainNode([2.0], [nothing, Message([3.0])])
        msg = ForneyLab.updateNodeMessage!(1, node, Array{Float64})
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => [1.5]
        # Forward message
        node = initializeFixedGainNode([2.0], [Message([3.0]), nothing])
        msg = ForneyLab.updateNodeMessage!(2, node, Array{Float64})
        @fact node.interfaces[2].message => msg
        @fact node.interfaces[2].message.value => [6.0]
    end

    context("FixedGainNode should propagate a univariate GaussianDistribution") do
        # The following tests on the update rules correspond to nodes 3 and 4 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        context("Univariate GaussianDistribution with (m,V) parametrization") do
            # Backward message
            node = initializeFixedGainNode([2.0], [nothing, Message(GaussianDistribution(m=3.0, V=5.0))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [1.5]
            @fact node.interfaces[1].message.value.V => reshape([1.25], 1, 1)
            # Forward message
            node = initializeFixedGainNode([2.0], [Message(GaussianDistribution(m=3.0, V=5.0)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [6.0]
            @fact node.interfaces[2].message.value.V => reshape([20.0], 1, 1)
        end
        context("Univariate GaussianDistribution with (m,V=0) parametrization") do
            # Backward message
            node = initializeFixedGainNode([2.0], [nothing, Message(GaussianDistribution(m=3.0, V=0.0))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [1.5]
            @fact node.interfaces[1].message.value.V => reshape([0.0], 1, 1)
            # Forward message
            node = initializeFixedGainNode([2.0], [Message(GaussianDistribution(m=3.0, V=0.0)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [6.0]
            @fact node.interfaces[2].message.value.V => reshape([0.0], 1, 1)
        end
        context("Univariate GaussianDistribution with (m,W) parametrization") do
            # Backward message
            node = initializeFixedGainNode([2.0], [nothing, Message(GaussianDistribution(m=3.0, W=2.0))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [1.5]
            @fact node.interfaces[1].message.value.W => reshape([8.0], 1, 1)
            # Forward message
            node = initializeFixedGainNode([2.0], [Message(GaussianDistribution(m=3.0, W=2.0)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [6.0]
            @fact node.interfaces[2].message.value.W => reshape([0.5], 1, 1)
        end
        context("Univariate GaussianDistribution with (xi,W) parametrization") do
            # Backward message
            node = initializeFixedGainNode([2.0], [nothing, Message(GaussianDistribution(xi=6.0, W=2.0))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.xi => [12.0]
            @fact node.interfaces[1].message.value.W => reshape([8.0], 1, 1)
            # Forward message
            node = initializeFixedGainNode([2.0], [Message(GaussianDistribution(xi=6.0, W=2.0)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.xi => [3.0]
            @fact node.interfaces[2].message.value.W => reshape([0.5], 1, 1)
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
            node = initializeFixedGainNode(A, [nothing, Message(GaussianDistribution(m=mean, V=variance))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => inv(A) * mean
            @fact node.interfaces[1].message.value.V => inv(A) * variance * inv(A)'
            # Forward message
            node = initializeFixedGainNode(A, [Message(GaussianDistribution(m=mean, V=variance)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => A * mean
            @fact node.interfaces[2].message.value.V => A * variance * A'
        end
        context("Multivariate GaussianDistribution with (m,V=0) parametrization") do
            mean = [1.0:3.0]
            variance = zeros(3, 3)
            # Backward message
            node = initializeFixedGainNode(A, [nothing, Message(GaussianDistribution(m=mean, V=variance))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => inv(A) * mean
            @fact node.interfaces[1].message.value.V => inv(A) * variance * inv(A)'
            # Forward message
            node = initializeFixedGainNode(A, [Message(GaussianDistribution(m=mean, V=variance)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => A * mean
            @fact node.interfaces[2].message.value.V => A * variance * A'
        end
        context("Multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = inv(reshape([   4.0, 3.0, 2.0,
                                        3.0, 4.0, 3.0,
                                        2.0, 3.0, 4.0], 3, 3))
            # Backward message
            node = initializeFixedGainNode(A, [nothing, Message(GaussianDistribution(m=mean, W=precision))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => inv(A) * mean
            @fact node.interfaces[1].message.value.W => A' * precision * A
            # Forward message
            node = initializeFixedGainNode(A, [Message(GaussianDistribution(m=mean, W=precision)), nothing])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => A * mean
            @fact node.interfaces[2].message.value.W => inv(A)' * precision * inv(A)
        end
    end
end