#####################
# Unit tests
#####################

facts("Clear message tests") do
    context("clearMessage!() and clearMessages!() should clear messages") do
        initializeAdditionNode()
        n(:add_node).i[:in1].message = Message(Gaussian())
        n(:add_node).i[:in2].message = Message(Gaussian())
        n(:add_node).i[:out].message = Message(Gaussian())

        clearMessage!(n(:add_node).i[:in2])
        @fact n(:add_node).i[:in2].message --> nothing
        @fact n(:add_node).i[:in1].message --> Message(Gaussian())
        @fact n(:add_node).i[:out].message --> Message(Gaussian())
        clearMessages!(n(:add_node))
        @fact n(:add_node).i[:in1].message --> nothing
        @fact n(:add_node).i[:out].message --> nothing
    end
end
