module TestForneyLab

using FactCheck
using ForneyLab

context("Supported message types") do
	facts("It should initialize a Gaussian message") do
		@fact GaussianMessage().V => [1.0]
		@fact GaussianMessage().m => [0.0]
	end

	facts("It should initiatize a multiplication parameter as message") do
		@fact GeneralMessage(1.0).value => 1.0
		@fact GeneralMessage([1.0, 2.0]).value => [1.0, 2.0]
	end
end

context("Supported node types") do
	context("Multiplication node") do
		facts("It should initialize a multiplication node") do
			@fact typeof(MultiplicationNode()) => MultiplicationNode
			@fact length(MultiplicationNode().interfaces) => 3
		end

		facts("It should propagate a message") do
			#TODO: implement calculations
		end
	end

	context("Constant node") do
		facts("It should initialize a constant node") do
			@fact typeof(ConstantNode()) => ConstantNode
			@fact length(ConstantNode().interfaces) => 1
		end

		facts("It should propagate a message") do
			node = ConstantNode(GaussianMessage())
			@fact node.interfaces[1].message => nothing
			calculatemessage!(1,node,Array(GaussianMessage,0),GaussianMessage())
			@fact typeof(node.interfaces[1].message) => GaussianMessage
		end
	end
end

context("Interface") do
	facts("It should carry a message") do
		node = ConstantNode()
		interface = Interface(node, GeneralMessage())
		@fact interface.message.value => 1.0

		# TODO: Creating an interface on a node does not yet propagate
		# the information about that interface to the node itself.
		@fact node.interface.message.value => 1.0
		@fact node.interfaces[1].message.value => 1.0
	end
	
	facts("It should connect to a partner") do
		node1 = MultiplicationNode()
		node1.multiplier.message = GeneralMessage()
		node1.source.message = GaussianMessage()
		# node1.interfaces[1].message = GeneralMessage() 
		# node1.interfaces[2].message = GaussianMessage() 

		node2 = ConstantNode()
		node2.interface.message = GeneralMessage()
		# node2.interfaces[1].message = GeneralMessage()

		# Connect constant to multiplier input
		node1.multiplier.partner = node2.interface
		node2.interface.partner = node1.multiplier
		# node1.interfaces[1].partner = node2.interfaces[1]
		# node2.interfaces[1].partner = node1.interfaces[1]

		# TODO: This is quite some code to just connect two nodes.
		# Connecting nodes should be easier.
		@fact node1.multiplier.partner => node2.interface
		@fact node2.interface.partner => node1.interface
	end

	facts("It should throw an error when the partner interface carries a different message type") do

	end
end

context("interfaces function") do

end

# context("Edges") do
# 	facts("It should initialize interfaces") do
# 		# Initialize a multiplication node with a Gaussian source and scalar multiplier
# 		node = MultiplicationNode();
# 		node.multiplier = Interface(node, nothing, GeneralMessage());
# 		node.source = Interface(node, nothing, GaussianMessage());
# 		node.sink = Interface(node);

# 		@fact typeof(node.multiplier.message) => GeneralMessage{Float64}
# 		@fact typeof(node.source.message) => GaussianMessage
# 		@fact typeof(node.sink.message) => Nothing
# 	end

# 	facts("It should initiatize edges as a combination of two interfaces") do
# 	end
# end

context("calculatemessage function") do
	
end

end # module TestForneyLab