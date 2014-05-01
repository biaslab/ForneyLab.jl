facts("ConstantNode") do
    context("ConstantNode() should initialize a ConstantNode with 1 interface") do
        node = ConstantNode()
        @fact typeof(node) => ConstantNode
        @fact length(node.interfaces) => 1
        @fact node.interface => node.interfaces[1]
    end

    context("ConstantNode should propagate a GaussianMessage") do
        node = ConstantNode(GaussianMessage(m=[2.0], V=[4.0]))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node, Array(None, 0))
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => GaussianMessage
        @fact node.interfaces[1].message.m => [2.0]
        @fact node.interfaces[1].message.V => reshape([4.0], 1, 1)
    end

    context("ConstantNode should propagate a GeneralMessage") do
        node = ConstantNode(GeneralMessage([1.0, 2.0]))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node, Array(None, 0))
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => GeneralMessage
        @fact node.interfaces[1].message.value => [1.0, 2.0]
    end
end