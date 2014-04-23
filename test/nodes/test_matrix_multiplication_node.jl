context("General node properties") do
	facts("MatrixMultiplicationNode should initialize a matrix multiplication node") do
		@fact typeof(MatrixMultiplicationNode()) => MatrixMultiplicationNode
		@fact length(MatrixMultiplicationNode().interfaces) => 2
		@fact typeof(MatrixMultiplicationNode(A=[1.0]).A) => Array{Float64,2} # cast single value to matrix
	end
end

context("Passing general messages") do
	facts("MatrixMultiplication node should propagate a general message") do
		node = MatrixMultiplicationNode([2.0])
		message = GeneralMessage([3.0])
		@fact calculatemessage!(1, node, [message, message]).value => [1.5] # backward propagation
		@fact calculatemessage!(2, node, [message, message]).value => [6] # forward propagation
	end
end

context("Gaussian message passing") do
	facts("MatrixMultiplication node should propagate a Gaussian message") do
		node = MatrixMultiplicationNode([2.0])
		message = GaussianMessage(m=[3.0], V=[5.0])
		backward_message = calculatemessage!(1, node, [message, message]) # backward propagation
		@fact backward_message.m => [1.5]
		@fact backward_message.V => reshape([1.25],1,1)
		forward_message = calculatemessage!(2, node, [message, message]) # forward propagation
		@fact forward_message.m => [6.0]
		@fact forward_message.V => reshape([50.0],1,1)
	end
end