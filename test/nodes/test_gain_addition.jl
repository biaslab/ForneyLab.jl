#####################
# Unit tests
#####################

facts("GainAdditionNode unit tests") do
    context("GainAdditionNode() should initialize a GainAdditionNode with 3 interfaces") do
        FactorGraph()
        GainAdditionNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:in1] --> n(:node).interfaces[1]
        @fact n(:node).i[:in2] --> n(:node).interfaces[2]
        @fact n(:node).i[:out] --> n(:node).interfaces[3]
        @fact typeof(n(:node).A) --> Array{Float64, 2}
    end

    msg_internal = Message(GaussianDistribution())
    context("GainAdditionNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result") do
        A = [2.0]

        # Forward
        initializeGainAdditionNode(A, [Message(GaussianDistribution(m=0.0, V=1.0)), Message(GaussianDistribution(m=1.0, V=2.0)), nothing])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(GaussianDistribution(m=0.0, V=1.0)), Message(GaussianDistribution(m=1.0, V=2.0)), nothing],
                                msg_internal.payload)

        # Backward
        initializeGainAdditionNode(A,[nothing, Message(GaussianDistribution(m=0.0, V=2.0)), Message(GaussianDistribution(m=1.0, V=2.0))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[1]))
        @fact msg_internal.payload --> GaussianDistribution(m=0.5, V=1.0)

        initializeGainAdditionNode(A, [Message(GaussianDistribution(m=0.0, V=2.0)), nothing, Message(GaussianDistribution(m=1.0, V=2.0))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(GaussianDistribution(m=0.0, V=2.0)), nothing, Message(GaussianDistribution(m=1.0, V=2.0))],
                                msg_internal.payload)
    end

    msg_internal = Message(MvGaussianDistribution())
    context("GainAdditionNode should be able to pass MvGaussianDistributions: using shortcut rules or internal graph should yield same result (m,V) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        initializeGainAdditionNode(A,[nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[1]))
        @fact msg_internal.payload --> MvGaussianDistribution(m=[0.8, -0.2], V=[1.56 -1.44; -1.44 1.56])

        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: using shortcut rules or internal graph should yield same result (m,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        initializeGainAdditionNode(A,[nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[1]))
        @fact msg_internal.payload --> MvGaussianDistribution(m=[0.8, -0.2], W=[8.0+(2/3) 8.0; 8.0 8.0+(2/3)])

        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: using shortcut rules or internal graph should yield same result (xi,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        initializeGainAdditionNode(A,[nothing, Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[1]))
        @fact msg_internal.payload --> MvGaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end

    context("GainAdditionNode should be able to pass MvGaussianDistributions: using shortcut rules or internal graph should yield same result (different parametrizations)") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[3]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                3,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing],
                                msg_internal.payload)

        # Backward
        initializeGainAdditionNode(A,[nothing, Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[1]))
        @fact msg_internal.payload --> MvGaussianDistribution(m=[0.4, -0.1], V=[0.78 -0.72; -0.72 0.78])

        initializeGainAdditionNode(A, [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = execute(ForneyLab.generateSumProductSchedule(n(:gac_node).interfaces[2]))
        FactorGraph()
        validateOutboundMessage(GainAdditionNode(A),
                                2,
                                [Message(MvGaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(MvGaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))],
                                msg_internal.payload)
    end
end

#####################
# Integration tests
#####################
