facts("SumProduct.collectInbounds() tests") do
    context("collectInbounds() should add the proper message/marginal") do
        # Standard
        (node, edges) = initializeGaussianNode()
        @fact SumProduct.collectInbounds(node.i[:out]) => (3, [node.i[:mean].partner.message, node.i[:precision].partner.message, nothing])

        # Composite node
        node = initializeGainEqualityNode(eye(1), Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        @fact SumProduct.collectInbounds(node.i[:out]) => (3, [node.i[:in1].partner.message, node.i[:in2].partner.message, nothing])
    end
end

# Test SumProduct specific functionality
include("test_generate_schedule.jl")
