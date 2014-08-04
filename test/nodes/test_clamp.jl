#####################
# Unit tests
#####################

facts("ClampNode unit tests") do
    context("ClampNode() should initialize a ClampNode with an out interface") do
        node = ForneyLab.ClampNode()
        @fact typeof(node) => ClampNode
        @fact typeof(node.out) => Interface
    end

    context("ClampNode should write a distribution to its interface") do
        payload = GaussianDistribution(m=2.0, V=4.0)
        node = ForneyLab.ClampNode(payload)
        @fact node.out.message => payload
    end

    context("ClampNode should write a message to its interface") do
        payload = Message(GaussianDistribution(m=2.0, V=4.0))
        node = ForneyLab.ClampNode(payload)
        @fact node.out.message => payload
    end
end
