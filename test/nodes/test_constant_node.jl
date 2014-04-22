context("General node properties") do
	facts("ConstantNode should initialize a constant node") do
		@fact typeof(ConstantNode()) => ConstantNode
		@fact length(ConstantNode().interfaces) => 1
	end
end

context("Sending out a message") do
	facts("Constant node should propagate a message") do
		node = ConstantNode(GaussianMessage(m=[2.0], V=[4.0]))
		@fact node.interfaces[1].message => nothing
		calculatemessage!(1,node,Array(GaussianMessage,0))
		@fact typeof(node.interfaces[1].message) => GaussianMessage
		@fact node.interfaces[1].message.m => [2.0]
		@fact node.interfaces[1].message.V => [4.0]
	end
end