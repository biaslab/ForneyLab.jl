#####################
# Unit tests
#####################

facts("Clear message tests") do
    context("clearMessage!() and clearMessages!() should clear messages") do
        initializeAdditionNode()
        n(:add_node).i[:in1].message = Message(GaussianDistribution())
        n(:add_node).i[:in2].message = Message(GaussianDistribution())
        n(:add_node).i[:out].message = Message(GaussianDistribution())

        clearMessage!(n(:add_node).i[:in2])
        @fact n(:add_node).i[:in2].message --> nothing
        @fact n(:add_node).i[:in1].message --> Message(GaussianDistribution())
        @fact n(:add_node).i[:out].message --> Message(GaussianDistribution())
        clearMessages!(n(:add_node))
        @fact n(:add_node).i[:in1].message --> nothing
        @fact n(:add_node).i[:out].message --> nothing
    end
end
