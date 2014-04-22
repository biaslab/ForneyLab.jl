# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl

module TestForneyLab

using FactCheck
using ForneyLab

context("General node properties") do
	facts("Node properties should include interfaces and name") do
		for(NodeType in subtypes(Node))
			@fact typeof(NodeType().interfaces) => Array{Interface,1} # Check for interface array
			@fact length(NodeType().interfaces) >= 1 => true # Check length of interface array
			@fact typeof(NodeType().name) => ASCIIString
		end
	end

	facts("Node constructor should assign a name") do
		for(NodeType in subtypes(Node))
			my_node = NodeType(;name="my_name")
			@fact my_node.name => "my_name"
		end
	end

	facts("Nodes should couple interfaces to themselves") do
		for(NodeType in subtypes(Node))
			my_node = NodeType()
			for(interface_id in 1:length(my_node.interfaces))
				# Check if the node interfaces couple back to the same node 
				@fact typeof(my_node.interfaces[interface_id].node) => typeof(my_node)
			end
		end
	end

	facts("Nodes should be able to calculate a message") do
		for(NodeType in subtypes(Node))
			# Check if method description contains node type
			@fact contains(string(methods(calculatemessage!)),string("::", NodeType)) => true
		end
	end
end

# Node and message specific tests are in separate files
include("test_messages.jl")
include("nodes/test_constant_node.jl")
include("nodes/test_multiplication_node.jl")
include("nodes/test_matrix_multiplication_node.jl")

# Helper function for initializing a pair of nodes
function initializepairofnodes()
	# Initialize some nodes
	node1 = MultiplicationNode()
	node1.interfaces[1].message = GaussianMessage() 
	node1.interfaces[2].message = GeneralMessage(1.0) 
	node2 = ConstantNode()
	node2.interfaces[1].message = GeneralMessage(2.0)
	return node1, node2
end

# Helper function for node comparison
function testinterfaceconnections(node1::MultiplicationNode, node2::ConstantNode)
	# Check that nodes are properly connected
	@fact node1.interfaces[2].message.value => 1.0
	@fact node2.interfaces[1].message.value => 2.0
	@fact node1.interfaces[2].partner.message.value => 2.0
	@fact node2.interfaces[1].partner.message.value => 1.0
	# Check that pointers are initiatized correctly
	@fact node1.in2.message.value => 1.0
	@fact node2.interface.message.value => 2.0
	@fact node1.in2.partner.message.value => 2.0
	@fact node2.interface.partner.message.value => 1.0
end

context("Connecting multiple nodes") do
	facts("Nodes can directly be coupled through interfaces by using the interfaces array") do
		(node1, node2) = initializepairofnodes()
		# Couple the interfaces that carry GeneralMessage
		node1.interfaces[2].partner = node2.interfaces[1]
		node2.interfaces[1].partner = node1.interfaces[2]
		testinterfaceconnections(node1, node2)
	end

	facts("Nodes can directly be coupled through interfaces by using the explicit interface names") do
		(node1, node2) = initializepairofnodes()
		# Couple the interfaces that carry GeneralMessage
		node1.in2.partner = node2.interface
		node2.interface.partner = node1.in2
		testinterfaceconnections(node1, node2)
	end

	facts("Nodes can be coupled by edges by using the interfaces array") do
		(node1, node2) = initializepairofnodes()
		# Couple the interfaces that carry GeneralMessage
		edge = Edge(node2.interfaces[1],node1.interfaces[2]) # Edge from node 2 to node 1
		testinterfaceconnections(node1, node2)
	end

	facts("Nodes can be coupled by edges using the explicit interface names") do
		(node1, node2) = initializepairofnodes()
		# Couple the interfaces that carry GeneralMessage
		edge = Edge(node2.interface,node1.in2) # Edge from node 2 to node 1
		testinterfaceconnections(node1, node2)
	end

	facts("Edge should throw an error when messages are of different types") do
		(node1, node2) = initializepairofnodes()
		# Couple the gaussian interface gaussian to the constant interface 
		@fact_throws Edge(node2.interfaces[1],node1.interfaces[1])
	end

	facts("Edge should throw an error when two interfaces on the same node are connected") do
		node = MultiplicationNode()
		node.interfaces[2].message = GaussianMessage() 
		node.interfaces[3].message = GaussianMessage()
		# Connect output directly to input 
		@fact_throws Edge(node.interfaces[3],node.interfaces[2])
	end

end

context("Message passing over interfaces") do
	facts("calculatemessage! should calculate an output message") do
		node1 = ConstantNode(GeneralMessage(3.0))
		node2 = MatrixMultiplicationNode([2.0])
		Edge(node1.interface, node2.in1)
		@fact node2.out.message => nothing
		# Request message on node for which the input is unknown
		calculatemessage!(node2.out)
		@fact typeof(node2.out.message) => GeneralMessage
		@fact node2.out.message.value => [6.0]
	end

	facts("calculatemessage! should recursively calculate a required message") do
		# Define three nodes in series
		node1 = ConstantNode(GeneralMessage(3.0))
		node2 = MatrixMultiplicationNode([2.0])
		node3 = MatrixMultiplicationNode([2.0])
		Edge(node1.interface, node2.in1)
		Edge(node2.out, node3.in1)
		@fact node3.out.message => nothing
		# Request message on node for which the input is unknown
		calculatemessage!(node3.out)
		@fact typeof(node3.out.message) => GeneralMessage
		@fact node3.out.message.value => [12.0]
	end

	facts("calculatemessage! should throw an error if input node and interface do not match") do
		(node1, node2) = initializepairofnodes()
		@fact_throws calculatemessage!(node1.out, node2)
	end

	facts("calculatemessage! should throw an error if the message-less interface has no partner") do
		node = MultiplicationNode()
		@fact_throws calculatemessage!(node.out)
	end
end

end # module TestForneyLab