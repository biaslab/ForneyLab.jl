facts("GainEqualityCompositeNode") do
    context("GainEqualityCompositeNode() should initialize a GainEqualityCompositeNode with 3 interfaces") do
        node = GainEqualityCompositeNode()
        @fact typeof(node) => GainEqualityCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.out1 => node.interfaces[2]
        @fact node.out2 => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainEqualityCompositeNode() should define an internal Equality and FixedGain node") do
        node = GainEqualityCompositeNode([5.0])
        @fact typeof(node.equality_node) => EqualityNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0],1,1)
    end

    context("GainEqualityCompositeNode() should point its own interfaces to the internal node interfaces") do
        node = GainEqualityCompositeNode([1.0])
        @fact node.interfaces[1] => node.equality_node.interfaces[1]
        @fact node.interfaces[2] => node.fixed_gain_node.interfaces[2]
        @fact node.interfaces[3] => node.equality_node.interfaces[3]
    end

    context("A GainEqualityCompositeNode should be able to pass a message through its internals") do
        #         _________
        #     in1 |       | out2
        # [N]-----|->[=]--|------>
        #         |   |   |
        #         |   v   |
        #         |  [A]  |
        #         |___|___|
        #             | out1
        #             v 
        #            [N]

        node = GainEqualityCompositeNode([2.0])
        c_node1 = ConstantNode(GaussianMessage(W=[1.0], xi=[1.0]))
        c_node2 = ConstantNode(GaussianMessage(W=[1.0], xi=[1.0]))
        Edge(c_node1.interfaces[1], node.interfaces[1])
        Edge(node.interfaces[2], c_node2.interfaces[1])
        msg = calculateMessage!(node.interfaces[3])
        println(msg)
    end
end