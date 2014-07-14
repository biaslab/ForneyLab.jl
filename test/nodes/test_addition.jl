#####################
# Unit tests
#####################

facts("AdditionNode unit tests") do
    context("AdditionNode() should initialize an AdditionNode with 3 interfaces") do
        node = AdditionNode()
        @fact typeof(node) => AdditionNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
    end
end

#####################
# Integration tests
#####################

facts("AdditionNode integration tests") do
    context("AdditionNode should add two Floats") do
        # Forward message
        node = initializeAdditionNode([Message(2.0), Message(3.0), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Float64)
        @fact node.interfaces[3].message => msg
        @fact node.interfaces[3].message.value => 5.0
        # Backward message
        node = initializeAdditionNode([nothing, Message(2.0), Message(3.0)])
        msg = ForneyLab.updateNodeMessage!(1, node, Float64)
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => 1.0
    end

    context("AdditionNode should add two Arrays") do
        # Forward message
        node = initializeAdditionNode([Message([1.0, 2.0]), Message([3.0, 4.0]), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, Array{Float64})
        @fact node.interfaces[3].message => msg
        @fact node.interfaces[3].message.value => [4.0, 6.0]
        # Backward message
        node = initializeAdditionNode([nothing, Message([1.0, 2.0]), Message([3.0, 4.0])])
        msg = ForneyLab.updateNodeMessage!(1, node, Array{Float64})
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => [2.0, 2.0]
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a univariate GaussianDistribution") do
        context("Univariate GaussianDistribution with (m,V) parametrization") do
            # Forward message
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], V=[2.0])), Message(GaussianDistribution(m=[3.0], V=[4.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.value.m => [4.0]
            @fact node.interfaces[3].message.value.V => reshape([6.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=[1.0], V=[2.0])), Message(GaussianDistribution(m=[3.0], V=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [2.0]
            @fact node.interfaces[1].message.value.V => reshape([6.0], 1, 1)
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], V=[2.0])), nothing, Message(GaussianDistribution(m=[3.0], V=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [2.0]
            @fact node.interfaces[2].message.value.V => reshape([6.0], 1, 1)
        end

        context("Univariate GaussianDistribution with (m,W) parametrization") do
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], W=[2.0])), Message(GaussianDistribution(m=[3.0], W=[4.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.value.m => [4.0]
            @fact node.interfaces[3].message.value.W => reshape([4.0/3.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=[1.0], W=[2.0])), Message(GaussianDistribution(m=[3.0], W=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [2.0]
            @fact node.interfaces[1].message.value.W => reshape([4.0/3.0], 1, 1)
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], W=[2.0])), nothing, Message(GaussianDistribution(m=[3.0], W=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [2.0]
            @fact node.interfaces[2].message.value.W => reshape([4.0/3.0], 1, 1)
        end

        context("Univariate GaussianDistribution with different parametrizations") do
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], V=[0.5])), Message(GaussianDistribution(m=[3.0], W=[4.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.value.m => [4.0]
            @fact node.interfaces[3].message.value.W => reshape([4.0/3.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=[1.0], V=[0.5])), Message(GaussianDistribution(m=[3.0], W=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [2.0]
            @fact node.interfaces[1].message.value.W => reshape([4.0/3.0], 1, 1)
            node = initializeAdditionNode([Message(GaussianDistribution(m=[1.0], V=[0.5])), nothing, Message(GaussianDistribution(m=[3.0], W=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [2.0]
            @fact node.interfaces[2].message.value.W => reshape([4.0/3.0], 1, 1)
        end

        context("Univariate GaussianDistribution with (xi,V) parametrization") do
            node = initializeAdditionNode([Message(GaussianDistribution(xi=[1.0], V=[2.0])), Message(GaussianDistribution(xi=[3.0], V=[4.0])), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.value.xi, [14/6]) => true
            @fact node.interfaces[3].message.value.V => reshape([6.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(xi=[1.0], V=[2.0])), Message(GaussianDistribution(xi=[3.0], V=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact isApproxEqual(node.interfaces[1].message.value.xi, [10/6]) => true
            @fact node.interfaces[1].message.value.V => reshape([6.0], 1, 1)
            node = initializeAdditionNode([Message(GaussianDistribution(xi=[1.0], V=[2.0])), nothing, Message(GaussianDistribution(xi=[3.0], V=[4.0]))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact isApproxEqual(node.interfaces[2].message.value.xi, [10/6]) => true
            @fact node.interfaces[2].message.value.V => reshape([6.0], 1, 1)
        end
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a multivariate GaussianDistribution") do
        context("Multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance)), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.value.m => [2.0, 4.0, 6.0]
            @fact node.interfaces[3].message.value.V => 2.0*variance
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [0.0, 0.0, 0.0]
            @fact node.interfaces[1].message.value.V => 2.0*variance
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, V=variance)), nothing, Message(GaussianDistribution(m=mean, V=variance))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [0.0, 0.0, 0.0]
            @fact node.interfaces[2].message.value.V => 2.0*variance
        end

        context("Multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision)), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.value.m, [2.0, 4.0, 6.0]) => true
            @fact isApproxEqual(node.interfaces[3].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.m => [0.0, 0.0, 0.0]
            @fact isApproxEqual(node.interfaces[1].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, W=precision)), nothing, Message(GaussianDistribution(m=mean, W=precision))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.m => [0.0, 0.0, 0.0]
            @fact isApproxEqual(node.interfaces[2].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
        end

        context("Multivariate GaussianDistribution with different parametrizations") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision))), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.value.m, [2.0, 4.0, 6.0]) => true
            @fact isApproxEqual(node.interfaces[3].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[1].message => msg
            @fact isApproxEqual(node.interfaces[1].message.value.m, [0.0, 0.0, 0.0]) => true
            @fact isApproxEqual(node.interfaces[1].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            node = initializeAdditionNode([Message(GaussianDistribution(m=mean, W=precision)), nothing, Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            ensureMWParametrization!(msg.value)
            @fact node.interfaces[2].message => msg
            @fact isApproxEqual(node.interfaces[2].message.value.m, [0.0, 0.0, 0.0]) => true
            @fact isApproxEqual(node.interfaces[2].message.value.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
        end

        context("Multivariate GaussianDistribution with (xi,V) parametrization") do
            xi = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance)), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianDistribution)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.value.xi, [1.0, 2.0, 3.0]) => true
            @fact node.interfaces[3].message.value.V => 2.0*variance
            # Backward messages
            node = initializeAdditionNode([nothing, Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianDistribution)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.value.xi => [0.0, 0.0, 0.0]
            @fact node.interfaces[1].message.value.V => 2.0*variance

            node = initializeAdditionNode([Message(GaussianDistribution(xi=xi, V=variance)), nothing, Message(GaussianDistribution(xi=xi, V=variance))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianDistribution)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.value.xi => [0.0, 0.0, 0.0]
            @fact node.interfaces[2].message.value.V => 2.0*variance
        end
    end
end