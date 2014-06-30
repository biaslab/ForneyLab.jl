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
            # TODO

            # Inverted
            node = initializeGaussianNode([GeneralMessage(2.0), InverseGammaMessage(a=3.0, b=1.0), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, InverseGammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the mean") do
            # Standard
            # TODO

            # Inverted
            node = initializeGaussianNode([nothing, InverseGammaMessage(a=3.0, b=1.0), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, InverseGammaMessage))
            @fact msg.m => [2.0]
            @fact msg.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the variance") do
            node = initializeGaussianNode([GeneralMessage(2.0), nothing, GeneralMessage(1.0)])
            msg = ForneyLab.updateNodeMessage!(2, node, GeneralMessage)
            @fact typeof(msg) => InverseGammaMessage
            @fact msg.a => -0.5
            @fact msg.b => 0.5
        end
    end

    context("Variational estimation") do
        context("Variational GaussianNode should propagate a forward message to y") do
        end

        context("Variational GaussianNode should propagate a backward variational message to the mean") do
            # Standard
            node = initializeVariationalGaussianNode([nothing, GammaMessage(a=3.0, b=1.0), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.W => reshape([3.0], 1, 1)
            # Inverse
            node = initializeVariationalGaussianNode([nothing, InverseGammaMessage(a=3.0, b=1.0), GeneralMessage(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, InverseGammaMessage))
            @fact typeof(msg) => GaussianMessage
            @fact msg.m => [2.0]
            @fact msg.V => reshape([4.0], 1, 1)
        end

        context("Variational GaussianNode should propagate a backward variational message to the variance or precision") do
            # Standard
            node = initializeVariationalGaussianNode([GaussianMessage(m=[4.0], W=[2.0]), nothing, GeneralMessage(2.0)])
            node.in2.edge.marginal = GammaMessage() # Let them know what I expect returned
            msg = ForneyLab.updateNodeMessage!(2, node, Union(GeneralMessage, GaussianMessage))
            @fact typeof(msg) => GammaMessage
            @fact msg.a => 1.5
            @fact msg.b => 2.25
            # Inverse
            node = initializeVariationalGaussianNode([GaussianMessage(m=[4.0], V=[1.0]), nothing, GeneralMessage(2.0)])
            node.in2.edge.marginal = InverseGammaMessage() # Let them know what I expect returned
            msg = ForneyLab.updateNodeMessage!(2, node, Union(GeneralMessage, GaussianMessage))
            @fact typeof(msg) => InverseGammaMessage
            @fact msg.a => -0.5
            @fact msg.b => 2.5
            # Throw error when marginal is not set and does not know what to return
            node = initializeVariationalGaussianNode([GaussianMessage(m=[4.0], V=[1.0]), nothing, GeneralMessage(2.0)])
            @fact_throws ForneyLab.updateNodeMessage!(2, node, Union(GeneralMessage, GaussianMessage))
        end
    end
end