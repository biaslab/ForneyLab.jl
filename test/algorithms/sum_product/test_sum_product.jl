facts("collectInbounds() tests") do
    context("collectInbounds() should collect the required inbound messages in an array") do
        # Standard
        initializeGaussianNode()
        @fact ForneyLab.collectInbounds(n(:node).i[:out], Val{symbol(sumProduct!)}) --> (3, [n(:node).i[:mean].partner.message, n(:node).i[:precision].partner.message, nothing])

        # Composite node
        initializeGainEqualityNode(eye(1), Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        @fact ForneyLab.collectInbounds(n(:gec_node).i[:out], Val{symbol(sumProduct!)}) --> (3, [n(:gec_node).i[:in1].partner.message, n(:gec_node).i[:in2].partner.message, nothing])
    end
end