facts("GainEqualityCompositeNode") do
    context("GainEqualityCompositeNode() should initialize a GainEqualityCompositeNode with 3 interfaces") do
        node = GainEqualityCompositeNode()
        @fact typeof(node) => GainEqualityCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.out1 => node.interfaces[2]
        @fact node.out2 => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2} # cast single value to matrix
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
end