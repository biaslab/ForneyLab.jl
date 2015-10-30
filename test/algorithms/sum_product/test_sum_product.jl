facts("SumProduct.collectInbounds() tests") do
    context("collectInbounds() should collect the required inbound messages in an array") do
        # Standard
        initializeGaussianNode()
        @fact SumProduct.collectInbounds(n(:node).i[:out]) --> (3, [n(:node).i[:mean].partner.message, n(:node).i[:precision].partner.message, nothing])


        # Include inbound message on outbound interface
        @fact SumProduct.collectInbounds(n(:node).i[:out], true) => (3, [n(:node).i[:mean].partner.message, n(:node).i[:precision].partner.message, n(:node).i[:out].partner.message])

        # Composite node
        initializeGainEqualityNode(eye(1), Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        @fact SumProduct.collectInbounds(n(:gec_node).i[:out]) --> (3, [n(:gec_node).i[:in1].partner.message, n(:gec_node).i[:in2].partner.message, nothing])
    end
end

# Test SumProduct specific functionality
include("test_generate_schedule.jl")
