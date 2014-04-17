context("Calculate messages for specific node types") do
	context("Multiplication node") do
		facts("MultiplicationNode should initialize a multiplication node") do
			@fact typeof(MultiplicationNode()) => MultiplicationNode
			@fact length(MultiplicationNode().interfaces) => 3
		end

		facts("Multiplication node should propagate a message") do
			#TODO: implement calculations
		end
	end

	context("Constant node") do
		facts("ConstantNode should initialize a constant node") do
			@fact typeof(ConstantNode()) => ConstantNode
			@fact length(ConstantNode().interfaces) => 1
		end

		facts("Constant node should propagate a message") do
			node = ConstantNode(GaussianMessage())
			@fact node.interfaces[1].message => nothing
			calculatemessage!(1,node,Array(GaussianMessage,0),GaussianMessage())
			@fact typeof(node.interfaces[1].message) => GaussianMessage
		end
	end
end
