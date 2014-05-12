# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl

module TestForneyLab

using FactCheck
using ForneyLab

# Helper function to check approximate equality
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < epsilon

facts("Helper functions") do
    context("ensureMatrix should convert an array with one element to a matrix type") do
        @fact typeof(ForneyLab.ensureMatrix([1.0])) => Array{Float64, 2} # Cast 1D to 2D array
        @fact ForneyLab.ensureMatrix([1.0]) => reshape([1.0], 1, 1)
        @fact ForneyLab.ensureMatrix(eye(2)) => eye(2)
    end
end

facts("General node properties") do
    context("Node properties should include interfaces and name") do
        for NodeType in subtypes(Node)
            @fact typeof(NodeType().interfaces) => Array{Interface, 1} # Check for interface array
            @fact length(NodeType().interfaces) >= 1 => true # Check length of interface array
            @fact typeof(NodeType().name) => ASCIIString
        end
    end

    context("Node constructor should assign a name") do
        for NodeType in subtypes(Node)
            my_node = NodeType(;name="my_name")
            @fact my_node.name => "my_name"
        end
    end

    context("Nodes should couple interfaces to themselves") do
        for NodeType in subtypes(Node)
            my_node = NodeType()
            for interface_id in 1:length(my_node.interfaces)
                # Check if the node interfaces couple back to the same node
                @fact my_node.interfaces[interface_id].node => my_node
            end
        end
    end

    context("Every node type should have at least 1 updateNodeMessage!() method") do
        for NodeType in subtypes(Node)
            # Check if method description contains node type
            @fact contains(string(methods(ForneyLab.updateNodeMessage!)), string("::", NodeType)) => true
        end
    end
end

# Node and message specific tests are in separate files
include("test_messages.jl")
include("nodes/test_addition.jl")
include("nodes/test_constant.jl")
include("nodes/test_equality.jl")
include("nodes/test_fixed_gain.jl")
include("nodes/composite/test_gain_addition.jl")
include("nodes/composite/test_gain_equality.jl")

# Helper function for initializing a pair of nodes
function initializePairOfNodes()
    # Initialize some nodes
    node1 = FixedGainNode()
    node1.interfaces[1].message = GeneralMessage(2.0) # Values differ to distinguish messages
    node1.interfaces[2].message = GeneralMessage(3.0)
    node2 = ConstantNode()
    node2.interfaces[1].message = GeneralMessage(1.0)
    return node1, node2
end

# Helper function for node comparison
function testInterfaceConnections(node1::FixedGainNode, node2::ConstantNode)
    # Check that nodes are properly connected
    @fact node1.interfaces[1].message.value => 2.0
    @fact node2.interfaces[1].message.value => 1.0
    @fact node1.interfaces[1].partner.message.value => 1.0
    @fact node2.interfaces[1].partner.message.value => 2.0
    # Check that pointers are initiatized correctly
    @fact node1.out.message.value => 3.0
    @fact node2.out.message.value => 1.0
    @fact node1.in1.partner.message.value => 1.0
    @fact node2.out.partner.message.value => 2.0
end

facts("Connections between nodes") do
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
        node = FixedGainNode()
        node.interfaces[1].message = GaussianMessage()
        node.interfaces[2].message = GaussianMessage()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

end

facts("Message passing over interfaces") do
    context("calculateMessage!() should return and write back an output message") do
        node1 = ConstantNode(GeneralMessage(3.0))
        node2 = FixedGainNode([2.0])
        Edge(node1.out, node2.in1)
        @fact node2.out.message => nothing
        # Request message on node for which the input is unknown
        msg = calculateMessage!(node2.out)
        @fact msg => node2.out.message # Returned message should be identical to message stored on interface.
        @fact typeof(node2.out.message) => GeneralMessage
        @fact node2.out.message.value => reshape([6.0], 1, 1)
        @fact node1.out.message_valid => false # consumed messages should be invalidated
        @fact node2.out.message_valid => true # fresh messages should be valid
    end

    context("calculateMessage!() should recursively calculate required inbound message") do
        # Define three nodes in series
        node1 = ConstantNode(GeneralMessage(3.0))
        node2 = FixedGainNode([2.0])
        node3 = FixedGainNode([2.0])
        Edge(node1.out, node2.in1)
        Edge(node2.out, node3.in1)
        @fact node3.out.message => nothing
        # Request message on node for which the input is unknown
        calculateMessage!(node3.out)
        @fact typeof(node3.out.message) => GeneralMessage
        @fact node3.out.message.value => reshape([12.0], 1, 1)
        @fact node1.out.message_valid => false # consumed messages should be invalidated
        @fact node2.out.message_valid => false
        @fact node3.out.message_valid => true # fresh messages should be valid
    end

    context("calculateMessage!() should throw an error if the specified interface does not belong to the specified node") do
        (node1, node2) = initializePairOfNodes()
        @fact_throws calculateMessage!(node1.out, node2)
    end

    context("calculateMessage!() should throw an error if one or more interfaces have no partner") do
        node = FixedGainNode()
        @fact_throws calculateMessage!(node.out)
    end

    context("calculateMarginal(edge) should check for legal forward/backward messages") do
        @fact_throws calculateMarginal(Edge())
        @fact_throws calculateMarginal(Edge(ConstantNode(), ConstantNode()))
    end

    context("calculateMarginal(forward_msg, backward_msg) should check equality of message types") do
        @fact_throws calculateMarginal(GaussianMessage(), GeneralMessage())
    end

    context("calculateMarginal(edge) should give correct result") do
        edge = Edge(ConstantNode(GaussianMessage(m=[0.0], V=[1.0])),
                    ConstantNode(GaussianMessage(m=[0.0], V=[1.0])))
        calculateForwardMessage!(edge)
        calculateBackwardMessage!(edge)
        marginal_msg = calculateMarginal(edge)
        ensureMVParametrization!(marginal_msg)
        @fact marginal_msg.m => [0.0]
        @fact maximum(abs(marginal_msg.V - reshape([0.5], 1, 1))) < epsilon => true
    end

    context("calculateMarginal(forward_msg, backward_msg) should give correct result") do
        marginal_msg = calculateMarginal(
                                GaussianMessage(m=[0.0], V=[1.0]),
                                GaussianMessage(m=[0.0], V=[1.0]))
        ensureMVParametrization!(marginal_msg)
        @fact marginal_msg.m => [0.0]
        @fact maximum(abs(marginal_msg.V - reshape([0.5], 1, 1))) < epsilon => true
    end
end

end # module TestForneyLab