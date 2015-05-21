#####################
# Unit tests
#####################

facts("GainAdditionNode unit tests") do
    context("GainAdditionNode() should initialize a GainAdditionNode with 3 interfaces") do
        FactorGraph()
        node = GainAdditionNode()
        @fact typeof(node) => GainAdditionNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    msg_internal = Message(GaussianDistribution())
    context("GainAdditionNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,V) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionNode(A,[nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.8, -0.2], V=[1.56 -1.44; -1.44 1.56])

        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionNode(A,[nothing, Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.8, -0.2], W=[8.0+(2/3) 8.0; 8.0 8.0+(2/3)])

        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (xi,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                3, 
                                [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionNode(A,[nothing, Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                2, 
                                [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (different parametrizations)") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionNode(A,[nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        node = initializeGainAdditionNode(A, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end
end

#####################
# Integration tests
#####################
