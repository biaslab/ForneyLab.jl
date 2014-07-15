#####################
# Unit tests
#####################

facts("ConstantNode unit tests") do
    context("ConstantNode() should initialize a ConstantNode with 1 interface") do
        node = ConstantNode()
        @fact typeof(node) => ConstantNode
        @fact length(node.interfaces) => 1
        @fact node.out => node.interfaces[1]
    end

    context("ConstantNode should throw an error when its value is set to a message") do
        @fact ConstantNode(1.0).value => 1.0
        @fact_throws ConstantNode(Message(1.0))
    end

    context("ConstantNode should propagate a GaussianMessage") do
        node = ConstantNode(GaussianDistribution(m=2.0, V=4.0))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => Message{GaussianDistribution}
        @fact node.interfaces[1].message.value.m => [2.0]
        @fact node.interfaces[1].message.value.V => reshape([4.0], 1, 1)
    end

    context("ConstantNode should propagate an arbitrary value") do
        node = ConstantNode([1.0, 2.0])
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
