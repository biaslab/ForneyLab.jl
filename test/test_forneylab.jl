# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl

module TestForneyLab

using FactCheck
using ForneyLab

include("test_style.jl") # Test style conventions on source files
include("integration_helpers.jl") # Helper file for integration tests, contains backgrounds and validations
include("test_helpers.jl") # Tests for ForneyLab helper methods

#####################
# Unit tests
#####################

facts("General node properties unit tests") do
    context("Node properties should include interfaces and name") do
        for NodeType in [subtypes(Node), subtypes(CompositeNode)]
            if NodeType != CompositeNode && NodeType != MockNode
                @fact typeof(NodeType().interfaces) => Array{Interface, 1} # Check for interface array
                @fact length(NodeType().interfaces) >= 1 => true # Check length of interface array
                @fact typeof(NodeType().name) => ASCIIString
            end
        end
    end

    context("Composite nodes should have property use_composite_update_rules") do
        for NodeType in [subtypes(CompositeNode)]
            @fact NodeType().use_composite_update_rules => true || false
        end
    end

    context("Node constructor should assign a name") do
        for NodeType in [subtypes(Node), subtypes(CompositeNode)]
            if NodeType != CompositeNode && NodeType != MockNode
                my_node = NodeType(;name="my_name")
                @fact my_node.name => "my_name"
            end
        end
    end

    context("Nodes should couple interfaces to themselves") do
        for NodeType in [subtypes(Node), subtypes(CompositeNode)]
            if NodeType != CompositeNode && NodeType != MockNode
                my_node = NodeType()
                for interface_id in 1:length(my_node.interfaces)
                    # Check if the node interfaces couple back to the same node
                    @fact my_node.interfaces[interface_id].node => my_node
                end
            end
        end
    end

    context("Every node type should have at least 1 updateNodeMessage!() method") do
        for NodeType in [subtypes(Node), subtypes(CompositeNode)]
            if NodeType != CompositeNode && NodeType != MockNode
                # Check if method description contains node type
                @fact contains(string(methods(ForneyLab.updateNodeMessage!)), string("::", NodeType)) => true
            end
        end
    end
end

facts("CalculateMessage!() unit tests") do
    context("calculateMessage!() should throw an error if the specified interface does not belong to the specified node") do
        (node1, node2) = initializePairOfNodes()
        @fact_throws calculateMessage!(node1.out, node2)
    end

    context("calculateMessage!() should throw an error if one or more interfaces have no partner") do
        node = FixedGainNode()
        @fact_throws calculateMessage!(node.out)
    end
end

facts("calculateMarginal unit tests") do
    context("setMarginal!() should preset a marginal") do
        (node1, node2) = initializePairOfMockNodes()
        edge = Edge(node1.out, node2.out)
        setMarginal!(edge, uninformative(GeneralMessage))
        @fact edge.head.message.value => 1.0
        @fact edge.tail.message.value => 1.0
        @fact edge.marginal.value => 1.0
    end

    context("calculateMarginal(forward_msg, backward_msg) should check equality of message types") do
        @fact_throws calculateMarginal(GaussianMessage(), GeneralMessage())
    end

    context("calculateMarginal!(edge) should give correct result and save the marginal to the edge") do
        edge = Edge(ConstantNode(GaussianMessage(m=[0.0], V=[1.0])),
                    ConstantNode(GaussianMessage(m=[0.0], V=[1.0])))
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
        marginal_msg = calculateMarginal!(edge)
        @fact edge.marginal => marginal_msg
        ensureMVParametrization!(marginal_msg)
        @fact marginal_msg.m => [0.0]
        @fact isApproxEqual(marginal_msg.V, reshape([0.5], 1, 1)) => true
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_msg = calculateMarginal(
                                GaussianMessage(m=[0.0], V=[1.0]),
                                GaussianMessage(m=[0.0], V=[1.0]))
        ensureMVParametrization!(marginal_msg)
        @fact marginal_msg.m => [0.0]
        @fact isApproxEqual(marginal_msg.V, reshape([0.5], 1, 1)) => true
    end
end

# Node and message specific tests are in separate files
include("test_messages.jl")
include("nodes/test_addition.jl")
include("nodes/test_constant.jl")
include("nodes/test_equality.jl")
include("nodes/test_fixed_gain.jl")
include("nodes/test_gaussian.jl")
include("nodes/composite/test_gain_addition.jl")
include("nodes/composite/test_gain_equality.jl")

#####################
# Integration tests
#####################

facts("Connections between nodes integration tests") do
    context("Nodes can directly be coupled through interfaces by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.interfaces[1].partner = node2.interfaces[1]
        node2.interfaces[1].partner = node1.interfaces[1]
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can directly be coupled through interfaces by using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.in1.partner = node2.out
        node2.out.partner = node1.in1
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can be coupled by edges by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can be coupled by edges using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.out, node1.in1) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Edge should throw an error when messages are of different types") do
        (node1, node2) = initializePairOfNodes()
        node1.interfaces[1].message = GaussianMessage()
        # Couple the gaussian interface gaussian to the constant interface
        @fact_throws Edge(node2.interfaces[1], node1.interfaces[1])
    end

    context("Edge should throw an error when two interfaces on the same node are connected") do
        node = initializeFixedGainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end
    
    context("Edge construction should couple interfaces to edge") do
        (node1, node2) = initializePairOfMockNodes()
        @fact node1.out.edge => nothing
        @fact node2.out.edge => nothing
        edge = Edge(node1.out, node2.out)
        @fact node1.out.edge => edge
        @fact node2.out.edge => edge
    end
end

facts("getAllNodesInGraph integration tests") do
    context("getAllNodesInGraph() should return an array of all nodes in the graph") do
        nodes = initializeLoopyGraph()
        found_nodes = ForneyLab.getAllNodesInGraph(nodes[1])
        @fact length(found_nodes) => length(nodes)
        for node in nodes
            @fact node in found_nodes => true
        end
    end
end

facts("calculateMessage!() integration tests") do
    context("calculateMessage!() should return and write back an output message") do
        (gain, constant) = initializePairOfNodes(A=[2.0], msg_gain_1=nothing, msg_gain_2=nothing, msg_const=GeneralMessage(3.0))
        Edge(constant.out, gain.in1)
        @fact gain.out.message => nothing
        # Request message on node for which the input is unknown
        msg = calculateMessage!(gain.out)
        @fact msg => gain.out.message # Returned message should be identical to message stored on interface.
        @fact typeof(gain.out.message) => GeneralMessage
        @fact gain.out.message.value => reshape([6.0], 1, 1)
    end

    context("calculateMessage!() should recursively calculate required inbound message") do
        # Define three nodes in series
        (node1, node2, node3) = initializeChainOfNodes()
        @fact node3.out.message => nothing
        # Request message on node for which the input is unknown
        calculateMessage!(node3.out)
        @fact typeof(node3.out.message) => GeneralMessage
        @fact node3.out.message.value => reshape([12.0], 1, 1)
    end

    context("calculateMessage!() should throw an error when there is an unbroken loop") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=[1.0], noise_V=[0.1])
        @fact_throws calculateMessage!(driver.out)
    end
end

facts("generateSchedule() and executeSchedule() integration tests") do
    (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=[1.0], noise_V=[0.1])

    context("generateSchedule() should throw an error when there is an unbroken loop") do
        @fact_throws generateSchedule(driver.out)
    end

    # Initial message
    setMessage!(add.in1, GaussianMessage(m=[2.0], V=[0.5]))
    setMessage!(add.out, GaussianMessage())

    context("generateSchedule() should auto-generate a feasible schedule") do
        # Generate schedule automatically
        schedule = generateSchedule(add.in2) # Message towards noise factor
        # All (but just) required calculations should be in the schedule
        @fact inhibitor.out in schedule => true
        @fact driver.out    in schedule => true
        @fact inhibitor.in1 in schedule => true
        @fact driver.in1    in schedule => true
        @fact add.in2       in schedule => true
        @fact add.in1       in schedule => false
        @fact add.out       in schedule => false
        @fact noise.out     in schedule => false
        # Validate correct relative order in schedule
        @fact findfirst(schedule, inhibitor.out)    < findfirst(schedule, driver.out)   => true
        @fact findfirst(schedule, driver.out)       < findfirst(schedule, add.in2)      => true
        @fact findfirst(schedule, driver.in1)       < findfirst(schedule, inhibitor.in1)=> true
        @fact findfirst(schedule, inhibitor.in1)    < findfirst(schedule, add.in2)      => true
    end

    context("generateSchedule() should correctly complete a partial schedule") do
        # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
        schedule = generateSchedule([driver.out, add.in2]) # Message towards noise factor
        # All (but just) required calculations should be in the schedule
        @fact schedule[1] => inhibitor.out
        @fact schedule[2] => driver.out
        @fact schedule[3] => driver.in1
        @fact schedule[4] => inhibitor.in1
        @fact schedule[5] => add.in2
    end

    context("executeSchedule() should correctly execute a schedule and return the result of the last step") do
        schedule = generateSchedule(add.in2)
        msg = ensureMVParametrization!(executeSchedule(schedule))
        @fact msg => add.in2.message
        @fact isApproxEqual(msg.m, [2.0]) => true
        @fact isApproxEqual(msg.V, reshape([1.5], 1, 1)) => true
    end

    context("executeSchedule() should accept edges") do
        (node1, node2, node3) = initializeChainOfNodes()
        schedule = [node1.out.edge, node2.out.edge]
        node1.out.message = GaussianMessage(W=[1.0], xi=[1.0]) 
        node1.out.partner.message = GaussianMessage(W=[1.0], xi=[1.0]) 
        node2.out.message = GaussianMessage(W=[1.0], xi=[1.0]) 
        node2.out.partner.message = GaussianMessage(W=[1.0], xi=[1.0])
        executeSchedule(schedule)
        @fact node1.out.edge.marginal.W => reshape([2.0], 1, 1)
        @fact node2.out.edge.marginal.W => reshape([2.0], 1, 1)
        @fact node1.out.edge.marginal.xi => [2.0]
        @fact node2.out.edge.marginal.xi => [2.0]
    end

    context("executeSchedule() should work as expeced in loopy graphs") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=[1.0], noise_V=[0.1])
        setMessage!(driver.out, GaussianMessage())
        schedule = generateSchedule(driver.out)
        for count = 1:100
            executeSchedule(schedule)
        end
        @fact typeof(driver.out.message) => GaussianMessage
        @fact ensureMVParametrization!(driver.out.message).m => [100.0] # For stop conditions at 100 cycles deep
    end

    context("executeSchedule() should be called repeatedly until convergence") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[1.1], B=[0.1], noise_m=[0.0], noise_V=[0.1])
        # Now set a breaker message and check that it works
        breaker_message = GaussianMessage(m=[10.0], V=[100.0])
        setMessage!(driver.out, breaker_message)
        prev_msg = deepcopy(breaker_message)
        converged = false
        schedule = generateSchedule(driver.out)
        while !converged
            msg = ensureMVParametrization!(executeSchedule(schedule))
            converged = isApproxEqual(prev_msg.m, msg.m)
            prev_msg = deepcopy(msg)
        end
        @fact isApproxEqual(driver.out.message.m, [0.0]) => true
    end
end

facts("clearMessage!(), clearMessages!(), clearAllMessages!() integration tests") do
    (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=[1.0], noise_V=[0.1])
    setMessage!(add.in1, GaussianMessage(m=[2.0], V=[0.5]))
    setMessage!(add.out, GaussianMessage())
    schedule = generateSchedule(add.in2)
    executeSchedule(schedule)
    clearMessage!(add.in2)
    @fact add.in2.message => nothing
    clearMessages!(add)
    @fact add.in1.message => nothing
    @fact add.out.message => nothing
    clearAllMessages!(add)
    for node in (driver, inhibitor, noise, add)
        for interface in node.interfaces
            @fact interface.message => nothing
        end
    end
end

try
    # Try to load user-defined extensions tests
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/test/test_forneylab_extensions.jl")
end

end # module TestForneyLab