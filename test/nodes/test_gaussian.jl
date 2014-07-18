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
            node = initializeGaussianNode([Message(2.0), Message(InverseGammaDistribution(a=3.0, b=1.0)), nothing])
            msg = ForneyLab.updateNodeMessage!(3, node, Union(Float64, InverseGammaDistribution), GaussianDistribution)
            @fact typeof(msg) => Message{GaussianDistribution}
            @fact msg.value.m => [2.0]
            @fact msg.value.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the mean") do
            # Standard
            # TODO

            # Inverted
            node = initializeGaussianNode([nothing, Message(InverseGammaDistribution(a=3.0, b=1.0)), Message(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(Float64, InverseGammaDistribution), GaussianDistribution)
            @fact msg.value.m => [2.0]
            @fact msg.value.V => reshape([0.5], 1, 1)
        end

        context("GaussianNode should propagate a backward message to the variance") do
            node = initializeGaussianNode([Message(2.0), nothing, Message(1.0)])
            msg = ForneyLab.updateNodeMessage!(2, node, Float64, InverseGammaDistribution)
            @fact typeof(msg) => Message{InverseGammaDistribution}
            @fact msg.value.a => -0.5
            @fact msg.value.b => 0.5
        end
    end


    context("Variational estimation") do
        context("Variational GaussianNode should propagate a forward message to y") do
        end

        context("Variational GaussianNode should propagate a backward variational message to the mean") do
            # Standard
            node = initializeVariationalGaussianNode([nothing, Message(GammaDistribution(a=3.0, b=1.0)), Message(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(Float64, GammaDistribution), GaussianDistribution)
            @fact typeof(msg) => Message{GaussianDistribution}
            @fact msg.value.m => [2.0]
            @fact msg.value.W => reshape([3.0], 1, 1)
            # Inverse
            node = initializeVariationalGaussianNode([nothing, Message(InverseGammaDistribution(a=3.0, b=1.0)), Message(2.0)])
            msg = ForneyLab.updateNodeMessage!(1, node, Union(Float64, InverseGammaDistribution), GaussianDistribution)
            @fact typeof(msg) => Message{GaussianDistribution}
            @fact msg.value.m => [2.0]
            @fact msg.value.V => reshape([4.0], 1, 1)
        end

        context("Variational GaussianNode should propagate a backward variational message to the variance or precision") do
            # Standard
            node = initializeVariationalGaussianNode([Message(GaussianDistribution(m=4.0, W=2.0)), nothing, Message(2.0)])
            node.in2.edge.marginal = GammaDistribution() # Let them know what I expect returned
            msg = ForneyLab.updateNodeMessage!(2, node, Union(Float64, GaussianDistribution), GammaDistribution)
            @fact typeof(msg) => Message{GammaDistribution}
            @fact msg.value.a => 1.5
            @fact msg.value.b => 2.25
            # Inverse
            node = initializeVariationalGaussianNode([Message(GaussianDistribution(m=4.0, V=1.0)), nothing, Message(2.0)])
            node.in2.edge.marginal = InverseGammaDistribution() # Let them know what I expect returned
            msg = ForneyLab.updateNodeMessage!(2, node, Union(Float64, GaussianDistribution), InverseGammaDistribution)
            @fact typeof(msg) => Message{InverseGammaDistribution}
            @fact msg.value.a => -0.5
            @fact msg.value.b => 2.5
        end
    end
end