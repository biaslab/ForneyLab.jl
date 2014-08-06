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

    context("AdditionNode() should allow parameters to be clamped") do
        # Fix in1
        node = AdditionNode(in1=GaussianDistribution())
        @fact typeof(node.in1.partner.node) => ForneyLab.ClampNode
        @fact node.in1.partner.message.value => GaussianDistribution()
        # Fix in2
        node = AdditionNode(in2=GaussianDistribution())
        @fact typeof(node.in2.partner.node) => ForneyLab.ClampNode
        @fact node.in2.partner.message.value => GaussianDistribution()
        # Fix out
        node = AdditionNode(out=GaussianDistribution())
        @fact typeof(node.out.partner.node) => ForneyLab.ClampNode
        @fact node.out.partner.message.value => GaussianDistribution()
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

    # TODO: fix these multivariate tests
    context("AdditionNode should propagate a multivariate GaussianDistribution") do
        context("AdditionNode should propagate a multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance)), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, V=variance)), nothing, Message(GaussianDistribution(m=mean, V=variance))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("AdditionNode should propagate a multivariate GaussianDistribution with (m,W) parametrization") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision)), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), nothing, Message(GaussianDistribution(m=mean, W=precision))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end

        context("AdditionNode should propagate a multivariate GaussianDistribution with (xi,V) parametrization") do
            xi = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance)), nothing],
                                    GaussianDistribution(xi=[1.0, 2.0, 3.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance))],
                                    GaussianDistribution(xi=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(xi=xi, V=variance)), nothing, Message(GaussianDistribution(xi=xi, V=variance))],
                                    GaussianDistribution(xi=[0.0, 0.0, 0.0], V=2.0*variance))
        end

        context("AdditionNode should propagate a multivariate GaussianDistribution with different parametrizations") do
            mean = [1.0:3.0]
            precision = reshape([4.0, 3.0, 2.0,
                                 3.0, 4.0, 3.0,
                                 2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision))), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    GaussianDistribution, 
                                    [nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    GaussianDistribution, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), nothing, Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end
    end
end

#####################
# Integration tests
#####################