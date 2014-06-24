#####################
# Unit tests
#####################

facts("GaussianNode unit tests") do
    context("GaussianNode() should initialize a GaussianNode with 3 interfaces") do
        node = GaussianNode()
        @fact typeof(node) => GaussianNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact node.variational => false # default variational to false
    end
end

#####################
# Integration tests
#####################

facts("GaussianNode integration tests") do
    context("Point estimates of y and m, so no approximation is required.") do
        context("GaussianNode should propagate a forward message to y") do
            # Standard
            node = initializeGaussianNode([GeneralMessage(2.0), GammaMessage(a=3.0, b=1.0, inverted=false), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.W => reshape([2.0], 1, 1)
            # Inverted
            node = initializeGaussianNode([GeneralMessage(2.0), GammaMessage(a=3.0, b=1.0, inverted=true), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the mean") do
            # Standard
            node = initializeGaussianNode([nothing, GammaMessage(a=3.0, b=1.0, inverted=false), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
            @fact msg.m => [2.0]
            @fact msg.W => reshape([2.0], 1, 1)
            # Inverted
            node = initializeGaussianNode([nothing, GammaMessage(a=3.0, b=1.0, inverted=true), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
            @fact msg.m => [2.0]
            @fact msg.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the variance") do
            node = initializeGaussianNode([GeneralMessage(2.0), nothing, GeneralMessage(1.0)])
            msg = ForneyLab.updateNodeMessage!(2, node, GeneralMessage)
            @fact typeof(msg) => GammaMessage
            @fact msg.a => -0.5
            @fact msg.b => 0.5
            @fact msg.inverted => true
        end
    end

    context("Variational estimation") do
        context("Variational GaussianNode should propagate a forward message to y") do
        end

        context("Variational GaussianNode should propagate a backward variational message to the mean") do
            # Inverted
            node = initializeVariationalGaussianNode([nothing, GammaMessage(a=3.0, b=1.0, inverted=true), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.V => reshape([4.0], 1, 1)
        end

        context("Variational GaussianNode should propagate a backward variational message to the variance") do
            # Inverted
            node = initializeVariationalGaussianNode([GaussianMessage(m=[4.0], V=[1.0]), nothing, GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(2, node, Union(GeneralMessage, GaussianMessage))
            @fact typeof(msg) => GammaMessage
            @fact msg.a => -0.5
            @fact msg.b => 2.5
            @fact msg.inverted => true
        end
    end
end