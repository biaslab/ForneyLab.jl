#####################
# Unit tests
#####################

facts("GainAdditionCompositeNode unit tests") do
    context("GainAdditionCompositeNode() should initialize a GainAdditionCompositeNode with 3 interfaces") do
        node = GainAdditionCompositeNode()
        @fact typeof(node) => GainAdditionCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainAdditionCompositeNode() should define an internal AdditionNode and FixedGainNode") do
        node = GainAdditionCompositeNode([5.0], false)
        @fact typeof(node.addition_node) => AdditionNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainAdditionCompositeNode() should point its own interfaces to the internal node interfaces") do
        node = GainAdditionCompositeNode([1.0], false)
        @fact node.in1.child => node.fixed_gain_node.interfaces[1]
        @fact node.in2.child => node.addition_node.interfaces[2]
        @fact node.out.child => node.addition_node.out
    end
end

#####################
# Integration tests
#####################

facts("GainAdditionCompositeNode integration tests") do
    context("Edge can connect a normal node to a GainAdditionCompositeNode") do
        (c_node, node) = initializeConstantAndGainAddNode()
        Edge(node.out, c_node.in2)
        @fact node.out.partner => c_node.in2 # Set correct partners
        @fact c_node.in2.partner => node.out
        @fact c_node.addition_node.in2.partner => node.out 
        @fact c_node.in2.child => c_node.addition_node.in2 # Set child 
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,V) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing])
        msg_shortcut = ForneyLab.updateNodeMessage!(3, node_gain_addition_composite, GaussianDistribution)
        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2))), nothing])
        msg_internal = ForneyLab.updateNodeMessage!(3, node_gain_addition_internal, GaussianDistribution)

        # TODO: msg_internal returns nothing

        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true

        # Backward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(1, node_gain_addition_composite, GaussianDistribution)
        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(1, node_gain_addition_internal, GaussianDistribution)
        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(2, node_gain_addition_composite, GaussianDistribution)
        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], V=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(2, node_gain_addition_internal, GaussianDistribution)
        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (m,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_shortcut = ForneyLab.updateNodeMessage!(3, node_gain_addition_composite, GaussianDistribution)
        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = ForneyLab.updateNodeMessage!(3, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true

        # Backward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(1, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(1, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(2, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(m=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(2, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
    end

    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (xi,W) parametrization") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_shortcut = ForneyLab.updateNodeMessage!(3, node_gain_addition_composite, GaussianDistribution)
        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = ForneyLab.updateNodeMessage!(3, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true

        # Backward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(1, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [nothing, Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(1, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(2, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(xi=[0.0, 0.0], W=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(2, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
    end
    context("GainAdditionCompositeNode should be able to pass GaussianDistributions: using shortcut rules or internal graph should yield same result (different parametrizations)") do
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)

        # Forward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_shortcut = ForneyLab.updateNodeMessage!(3, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2))), nothing])
        msg_internal = ForneyLab.updateNodeMessage!(3, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true

        # Backward
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(1, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [nothing, Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(1, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
        #Shortcut rules
        node_gain_addition_composite = initializeGainAdditionCompositeNode(A, true, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_shortcut = ForneyLab.updateNodeMessage!(2, node_gain_addition_composite, GaussianDistribution)        #Internal graph
        node_gain_addition_internal = initializeGainAdditionCompositeNode(A, false, [Message(GaussianDistribution(m=[0.0, 0.0], V=eye(2,2))), nothing, Message(GaussianDistribution(xi=[1.0, 2.0], W=2.0*eye(2,2)))])
        msg_internal = ForneyLab.updateNodeMessage!(2, node_gain_addition_internal, GaussianDistribution)        # Messages must be equal
        ensureMVParametrization!(msg_shortcut.value)
        ensureMVParametrization!(msg_internal.value)
        @fact isApproxEqual(msg_internal.value.m, msg_shortcut.value.m) => true
        @fact isApproxEqual(msg_internal.value.V, msg_shortcut.value.V) => true
    end
end