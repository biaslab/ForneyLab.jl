#####################
# Unit tests
#####################

facts("AdditionNode unit tests") do
    context("AdditionNode() should initialize an AdditionNode with 3 interfaces") do
        node = AdditionNode()
        @fact typeof(node) => AdditionNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
    end

    context("AdditionNode should add two Floats") do
        # Forward message
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                Float64, 
                                [Message(2.0), Message(3.0), nothing],
                                5.0)
        # Backward message
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                Float64, 
                                [nothing, Message(2.0), Message(3.0)],
                                1.0)
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                Float64, 
                                [Message(2.0), nothing, Message(3.0)],
                                1.0)
    end

    context("AdditionNode should add two Arrays") do
        # Forward message
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                Vector{Float64}, 
                                [Message([1.0, 2.0]), Message([3.0, 4.0]), nothing],
                                [4.0, 6.0])
        # Backward message
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                Vector{Float64}, 
                                [nothing, Message([1.0, 2.0]), Message([3.0, 4.0])],
                                [2.0, 2.0])
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                Vector{Float64}, 
                                [Message([1.0, 2.0]), nothing, Message([3.0, 4.0])],
                                [2.0, 2.0])
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a univariate GaussianDistribution") do
        context("AdditionNode should propagate a univariate GaussianDistribution with (m,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0)), nothing],
                                    GaussianDistribution(m=4.0, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), nothing, Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))        
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with (m,W) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))  
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with (xi,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0)), nothing],
                                    GaussianDistribution(xi=14/6, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), nothing, Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))  
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with different parametrizations") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))  
        end
    end

    # TODO: add multivariate tests
end

#####################
# Integration tests
#####################