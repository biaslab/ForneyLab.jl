#####################
# Unit tests
#####################

facts("TerminalNode unit tests") do
    context("TerminalNode() should initialize a TerminalNode with 1 interface") do
        node = TerminalNode()
        @fact typeof(node) => TerminalNode
        @fact length(node.interfaces) => 1
        @fact node.out => node.interfaces[1]
    end

    context("TerminalNode should throw an error when its value is set to a message") do
        @fact TerminalNode(1.0).value => 1.0
        @fact_throws TerminalNode(Message(1.0))
    end

    context("TerminalNode should propagate a GaussianMessage") do
        node = TerminalNode(GaussianDistribution(m=2.0, V=4.0))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => Message{GaussianDistribution}
        @fact node.interfaces[1].message.value.m => [2.0]
        @fact node.interfaces[1].message.value.V => reshape([4.0], 1, 1)
    end

    context("TerminalNode should propagate an arbitrary value") do
        node = TerminalNode([1.0, 2.0])
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => Message{Array{Float64, 1}}
        @fact node.interfaces[1].message.value => [1.0, 2.0]
    end
end

#####################
# Integration tests
#####################
