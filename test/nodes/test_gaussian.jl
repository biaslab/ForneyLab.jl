#####################
# Unit tests
#####################

facts("GaussianNode unit tests") do
    context("GaussianNode() should initialize a GaussianNode with 3 interfaces") do
        node = GaussianNode()
        @fact typeof(node) => GaussianNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
    end
end

#####################
# Integration tests
#####################

facts("GaussianNode integration tests") do
	context("Point estimates of y and m, so no approximation is required.") do
	    context("GaussianNode should propagate a forward message to y") do
	    	# Standard
	    	node = initializeGaussianNode([GeneralMessage(), GammaMessage(a=1.0, b=1.0, inverted=false), nothing])
  	        msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
  	        @fact true => false
  	        # Inverted
	    	node = initializeGaussianNode([GeneralMessage(), GammaMessage(a=1.0, b=1.0, inverted=true), nothing])
  	        msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
  	        @fact true => false
	    end

	    context("GaussianNode should propagate a backward message to the mean") do
	    	# Standard
	    	node = initializeGaussianNode([nothing, GammaMessage(a=1.0, b=1.0, inverted=false), GeneralMessage()])
  	        msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
  	        @fact true => false
  	        # Inverted
	    	node = initializeGaussianNode([nothing, GammaMessage(a=1.0, b=1.0, inverted=true), GeneralMessage()])
  	        msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
  	        @fact true => false
	    end

	    context("GaussianNode should propagate a backward message to the variance") do
	    	node = initializeGaussianNode([GeneralMessage(2.0), nothing, GeneralMessage(1.0)])
  	        msg = ForneyLab.updateNodeMessage!(2, node, GeneralMessage)
  	        @fact true => false
	    end
	end
end