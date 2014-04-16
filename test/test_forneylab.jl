module TestForneyLab

using FactCheck
using ForneyLab

context("Messages can be of different types") do
	facts("It should initialize a Gaussian message") do
		@fact ForneyLab.GaussianMessage().V => [1.0]
		@fact ForneyLab.GaussianMessage().m => [0.0]
	end

	facts("It should initiatize a multiplication parameter as message") do
		@fact ForneyLab.ScalarParameterMessage(1.0).value => 1.0
		@fact ForneyLab.ScalarParameterMessage([1.0, 2.0]).value => [1.0, 2.0]
	end
end

context("Nodes can be of different types") do
	facts("It should initialize a multiplication node") do
		@fact typeof(ForneyLab.MultiplicationNode()) => ForneyLab.MultiplicationNode
	end

	facts("It should initialize a constant node") do
		@fact typeof(ForneyLab.ConstantNode(2)) => ForneyLab.ConstantNode
	end
end

context("Interface") do
	facts("It should carry a message") do
		node = ForneyLab.MultiplicationNode()
		@fact ForneyLab.Interface(node, ForneyLab.GaussianMessage()).message.V => 1.0 
	end
	
	facts("It should connect to a partner") do
		node1 = ForneyLab.MultiplicationNode()
		node1.multiplier.message = ForneyLab.ScalarParameterMessage()
		node1.source.message = ForneyLab.GaussianMessage()

		node2 = ForneyLab.ConstantNode()
		node2.constant.message = ScalarParameterMessage()

		# TODO: build a connect function?
		node1.multiplier.partner = node2.interface
		node2.interface.partner = node1.multiplier
	end

	facts("It should throw an error when the partner interface carries a different message type") do
		node1 = ForneyLab.MultiplicationNode()
		node2 = ForneyLab.ConstantNode(2)
		# Connect nodes
	end
end

context("interfaces function") do

end

# context("Edges") do
# 	facts("It should initialize interfaces") do
# 		# Initialize a multiplication node with a Gaussian source and scalar multiplier
# 		node = ForneyLab.MultiplicationNode();
# 		node.multiplier = ForneyLab.Interface(node, nothing, ForneyLab.ScalarParameterMessage());
# 		node.source = ForneyLab.Interface(node, nothing, ForneyLab.GaussianMessage());
# 		node.sink = ForneyLab.Interface(node);

# 		@fact typeof(node.multiplier.message) => ForneyLab.ScalarParameterMessage{Float64}
# 		@fact typeof(node.source.message) => ForneyLab.GaussianMessage
# 		@fact typeof(node.sink.message) => Nothing
# 	end

# 	facts("It should initiatize edges as a combination of two interfaces") do
# 	end
# end

context("calculatemessage function") do
	
end

end # module TestForneyLab