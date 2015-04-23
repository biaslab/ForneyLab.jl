#####################
# Unit tests
#####################

facts("GainAdditionCompositeNode unit tests") do
    context("GainAdditionCompositeNode() should initialize a GainAdditionCompositeNode with 3 interfaces") do
        FactorGraph()
        node = GainAdditionCompositeNode()
        @fact typeof(node) => GainAdditionCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainAdditionCompositeNode() should define an internal AdditionNode and FixedGainNode") do
        FactorGraph()
        node = GainAdditionCompositeNode([5.0], false)
        @fact typeof(node.addition_node) => AdditionNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainAdditionCompositeNode() should point its own interfaces to the internal node interfaces") do
        FactorGraph()
        node = GainAdditionCompositeNode([1.0], false)
        @fact node.in1.child => node.fixed_gain_node.interfaces[1]
        @fact node.in2.child => node.addition_node.interfaces[2]
        @fact node.out.child => node.addition_node.out
    end

    msg_internal = Message(GaussianDistribution())
    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,V) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.8, -0.2], V=[1.56 -1.44; -1.44 1.56])

        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.8, -0.2], W=[8.0+(2/3) 8.0; 8.0 8.0+(2/3)])

        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (xi,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                3, 
                                [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                2, 
                                [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (different parametrizations)") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                3, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        node = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[1]))
        @fact msg_internal.payload => GaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        node = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(SumProduct.generateSchedule(node.interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionCompositeNode(A, true), 
                                2, 
                                [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end
end

#####################
# Integration tests
#####################

facts("GainAdditionCompositeNode integration tests") do
    context("Edge can connect a normal node to a GainAdditionCompositeNode") do
        (c_node, node) = initializeTerminalAndGainAddNode()
        Edge(node.out, c_node.in2)
        @fact node.out.partner => c_node.in2 # Set correct partners
        @fact c_node.in2.partner => node.out
        @fact c_node.addition_node.in2.partner => node.out 
        @fact c_node.in2.child => c_node.addition_node.in2 # Set child 
    end
end