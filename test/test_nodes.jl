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

	context("Matrix multiplication node") do
		facts("MatrixMultiplicationNode should initialize a matrix multiplication node") do
			@fact typeof(MatrixMultiplicationNode()) => MatrixMultiplicationNode
			@fact length(MatrixMultiplicationNode().interfaces) => 2
		end

		facts("MatrixMultiplication node should propagate a message") do
			#TODO: implement calculations
		end
	end

	context("Constant node") do
		facts("ConstantNode should initialize a constant node") do
			@fact typeof(ConstantNode()) => ConstantNode
			@fact length(ConstantNode().interfaces) => 1
		end

		facts("Constant node should propagate a message") do
			node = ConstantNode(GaussianMessage(m=[2.0], V=[4.0]))
			@fact node.interfaces[1].message => nothing
			calculatemessage!(1,node,Array(GaussianMessage,0))
			@fact typeof(node.interfaces[1].message) => GaussianMessage
			@fact node.interfaces[1].message.m => [2.0]
			@fact node.interfaces[1].message.V => [4.0]
		end
	end
end
