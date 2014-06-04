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
    context("AdditionNode should add two GeneralMessages") do
        # Forward message
        node = initializeAdditionNode([GeneralMessage(2.0), GeneralMessage(3.0), nothing])
        msg = ForneyLab.updateNodeMessage!(3, node, GeneralMessage)
        @fact node.interfaces[3].message => msg
        @fact node.interfaces[3].message.value => 5.0
        # Backward message
        node = initializeAdditionNode([nothing, GeneralMessage(2.0), GeneralMessage(3.0)])
        msg = ForneyLab.updateNodeMessage!(1, node, GeneralMessage)
        @fact node.interfaces[1].message => msg
        @fact node.interfaces[1].message.value => 1.0
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a univariate GaussianMessage") do
        context("Univariate GaussianMessage with (m,V) parametrization") do
            # Forward message
            node = initializeAdditionNode([GaussianMessage(m=[1.0], V=[2.0]), GaussianMessage(m=[3.0], V=[4.0]), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.m => [4.0]
            @fact node.interfaces[3].message.V => reshape([6.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=[1.0], V=[2.0]), GaussianMessage(m=[3.0], V=[4.0])])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [2.0]
            @fact node.interfaces[1].message.V => reshape([6.0], 1, 1)
            node = initializeAdditionNode([GaussianMessage(m=[1.0], V=[2.0]), nothing, GaussianMessage(m=[3.0], V=[4.0])])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [2.0]
            @fact node.interfaces[2].message.V => reshape([6.0], 1, 1)
        end

        context("Univariate GaussianMessage with (m,W) parametrization") do
            node = initializeAdditionNode([GaussianMessage(m=[1.0], W=[2.0]), GaussianMessage(m=[3.0], W=[4.0]), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.m => [4.0]
            @fact node.interfaces[3].message.W => reshape([4.0/3.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=[1.0], W=[2.0]), GaussianMessage(m=[3.0], W=[4.0])])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [2.0]
            @fact node.interfaces[1].message.W => reshape([4.0/3.0], 1, 1)
            node = initializeAdditionNode([GaussianMessage(m=[1.0], W=[2.0]), nothing, GaussianMessage(m=[3.0], W=[4.0])])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [2.0]
            @fact node.interfaces[2].message.W => reshape([4.0/3.0], 1, 1)
        end

        context("Univariate GaussianMessage with different parametrizations") do
            node = initializeAdditionNode([GaussianMessage(m=[1.0], V=[0.5]), GaussianMessage(m=[3.0], W=[4.0]), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.m => [4.0]
            @fact node.interfaces[3].message.W => reshape([4.0/3.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=[1.0], V=[0.5]), GaussianMessage(m=[3.0], W=[4.0])])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [2.0]
            @fact node.interfaces[1].message.W => reshape([4.0/3.0], 1, 1)
            node = initializeAdditionNode([GaussianMessage(m=[1.0], V=[0.5]), nothing, GaussianMessage(m=[3.0], W=[4.0])])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [2.0]
            @fact node.interfaces[2].message.W => reshape([4.0/3.0], 1, 1)
        end

        context("Univariate GaussianMessage with (xi,V) parametrization") do
            node = initializeAdditionNode([GaussianMessage(xi=[1.0], V=[2.0]), GaussianMessage(xi=[3.0], V=[4.0]), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.xi, [14/6]) => true
            @fact node.interfaces[3].message.V => reshape([6.0], 1, 1)
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(xi=[1.0], V=[2.0]), GaussianMessage(xi=[3.0], V=[4.0])])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact isApproxEqual(node.interfaces[1].message.xi, [10/6]) => true
            @fact node.interfaces[1].message.V => reshape([6.0], 1, 1)
            node = initializeAdditionNode([GaussianMessage(xi=[1.0], V=[2.0]), nothing, GaussianMessage(xi=[3.0], V=[4.0])])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact isApproxEqual(node.interfaces[2].message.xi, [10/6]) => true
            @fact node.interfaces[2].message.V => reshape([6.0], 1, 1)
        end
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a multivariate GaussianMessage") do
        context("Multivariate GaussianMessage with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([GaussianMessage(m=mean, V=variance), GaussianMessage(m=mean, V=variance), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact node.interfaces[3].message.m => [2.0, 4.0, 6.0]
            @fact node.interfaces[3].message.V => 2.0*variance
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=mean, V=variance), GaussianMessage(m=mean, V=variance)])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [0.0, 0.0, 0.0]
            @fact node.interfaces[1].message.V => 2.0*variance
            node = initializeAdditionNode([GaussianMessage(m=mean, V=variance), nothing, GaussianMessage(m=mean, V=variance)])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [0.0, 0.0, 0.0]
            @fact node.interfaces[2].message.V => 2.0*variance
        end

        context("Multivariate GaussianMessage with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([GaussianMessage(m=mean, W=precision), GaussianMessage(m=mean, W=precision), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.m, [2.0, 4.0, 6.0]) => true
            @fact isApproxEqual(node.interfaces[3].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=mean, W=precision), GaussianMessage(m=mean, W=precision)])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.m => [0.0, 0.0, 0.0]
            @fact isApproxEqual(node.interfaces[1].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            node = initializeAdditionNode([GaussianMessage(m=mean, W=precision), nothing, GaussianMessage(m=mean, W=precision)])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.m => [0.0, 0.0, 0.0]
            @fact isApproxEqual(node.interfaces[2].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
        end

        context("Multivariate GaussianMessage with different parametrizations") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([GaussianMessage(m=mean, W=precision), GaussianMessage(xi=precision*mean, V=inv(precision)), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.m, [2.0, 4.0, 6.0]) => true
            @fact isApproxEqual(node.interfaces[3].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(m=mean, W=precision), GaussianMessage(xi=precision*mean, V=inv(precision))])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[1].message => msg
            @fact isApproxEqual(node.interfaces[1].message.m, [0.0, 0.0, 0.0]) => true
            @fact isApproxEqual(node.interfaces[1].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
            node = initializeAdditionNode([GaussianMessage(m=mean, W=precision), nothing, GaussianMessage(xi=precision*mean, V=inv(precision))])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            ensureMWParametrization!(msg)
            @fact node.interfaces[2].message => msg
            @fact isApproxEqual(node.interfaces[2].message.m, [0.0, 0.0, 0.0]) => true
            @fact isApproxEqual(node.interfaces[2].message.W, reshape([2.0, 1.5, 1.0, 1.5, 2.0, 1.5, 1.0, 1.5, 2.0], 3, 3)) => true
        end

        context("Multivariate GaussianMessage with (xi,V) parametrization") do
            xi = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            node = initializeAdditionNode([GaussianMessage(xi=xi, V=variance), GaussianMessage(xi=xi, V=variance), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, GaussianMessage)
            @fact node.interfaces[3].message => msg
            @fact isApproxEqual(node.interfaces[3].message.xi, [1.0, 2.0, 3.0]) => true
            @fact node.interfaces[3].message.V => 2.0*variance
            # Backward messages
            node = initializeAdditionNode([nothing, GaussianMessage(xi=xi, V=variance), GaussianMessage(xi=xi, V=variance)])
            msg = ForneyLab.updateNodeMessage!(1, node, GaussianMessage)
            @fact node.interfaces[1].message => msg
            @fact node.interfaces[1].message.xi => [0.0, 0.0, 0.0]
            @fact node.interfaces[1].message.V => 2.0*variance

            node = initializeAdditionNode([GaussianMessage(xi=xi, V=variance), nothing, GaussianMessage(xi=xi, V=variance)])
            msg = ForneyLab.updateNodeMessage!(2, node, GaussianMessage)
            @fact node.interfaces[2].message => msg
            @fact node.interfaces[2].message.xi => [0.0, 0.0, 0.0]
            @fact node.interfaces[2].message.V => 2.0*variance
        end
    end
end