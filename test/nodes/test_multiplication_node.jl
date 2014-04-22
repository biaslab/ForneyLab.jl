context("General node properties") do
	facts("MultiplicationNode should initialize a multiplication node") do
		@fact typeof(MultiplicationNode()) => MultiplicationNode
		@fact length(MultiplicationNode().interfaces) => 3
	end
end

context("Passing general messages") do
	facts("Multiplication node should propagate a general message") do
		#TODO: implement calculations
	end
end

context("Gaussian message passing") do
	facts("Multiplication node should propagate a Gaussian message") do
		#TODO: implement calculations
	end
end

