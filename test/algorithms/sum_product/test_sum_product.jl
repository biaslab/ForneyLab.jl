facts("SumProduct.collectInbounds() tests") do
    context("collectInbounds() should collect the required inbound messages in an array") do
        # Standard
        (node, edges) = initializeGaussianNode()
        @fact SumProduct.collectInbounds(node.i[:out]) => (3, [node.i[:mean].partner.message, node.i[:precision].partner.message, nothing])

        # Include inbound message on outbound interface
        @fact SumProduct.collectInbounds(node.i[:out], true) => (3, [node.i[:mean].partner.message, node.i[:precision].partner.message, node.i[:out].partner.message])
    end
end

# Test SumProduct specific functionality
include("test_generate_schedule.jl")
