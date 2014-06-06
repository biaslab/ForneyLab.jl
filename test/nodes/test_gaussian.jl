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
	    	node = initializeGaussianNode([GeneralMessage(), GammaMessage(false), nothing])
  	        msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
  	        show(msg)
  	        # Inverted
	    	node = initializeGaussianNode([GeneralMessage(), GammaMessage(true), nothing])
  	        msg = ForneyLab.updateNodeMessage!(3, node, Union(GeneralMessage, GammaMessage))
  	        show(msg)
	    end

	    context("GaussianNode should propagate a backward message to the mean") do
	    	# Standard
	    	node = initializeGaussianNode([nothing, GammaMessage(false), GeneralMessage()])
  	        msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
  	        show(msg)
  	        # Inverted
	    	node = initializeGaussianNode([nothing, GammaMessage(true), GeneralMessage()])
  	        msg = ForneyLab.updateNodeMessage!(1, node, Union(GeneralMessage, GammaMessage))
  	        show(msg)
	    end

	    context("GaussianNode should propagate a backward message to the variance") do
	    	node = initializeGaussianNode([GeneralMessage(), nothing, GeneralMessage()])
  	        msg = ForneyLab.updateNodeMessage!(2, node, GeneralMessage)
  	        show(msg)
	    end
	end
end