facts("SumProduct.collectInbounds() tests") do
    context("collectInbounds() should add the proper message/marginal") do
        # Standard
        (node, edges) = initializeGaussianNode()
        @fact SumProduct.collectInbounds(node.out) => (3, [node.mean.partner.message, node.precision.partner.message, nothing])

        # Composite node
        node = initializeGainEqualityCompositeNode(eye(1), true, Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        @fact SumProduct.collectInbounds(node.out) => (3, [node.in1.partner.message, node.in2.partner.message, nothing])
    end
end

# Test SumProduct specific functionality
include("test_generate_schedule.jl")
