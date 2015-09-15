#####################
# Unit tests
#####################

facts("AdditionNode unit tests") do
    context("AdditionNode() should initialize an AdditionNode with 3 interfaces") do
        FactorGraph()
        AdditionNode(id=:node)
        @fact typeof(n(:node)) => AdditionNode
        @fact length(n(:node).interfaces) => 3
        @fact n(:node).i[:in1] => n(:node).interfaces[1]
        @fact n(:node).i[:in2] => n(:node).interfaces[2]
        @fact n(:node).i[:out] => n(:node).interfaces[3]
    end

    context("AdditionNode should add two DeltaDistribution{Float}") do
        # Forward message
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                [Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0)), nothing],
                                DeltaDistribution(5.0))
        # Backward message
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                [nothing, Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))],
                                DeltaDistribution(1.0))
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                [Message(DeltaDistribution(2.0)), nothing, Message(DeltaDistribution(3.0))],
                                DeltaDistribution(1.0))
    end

    context("AdditionNode should add two DeltaDistribution{Array}") do
        # Forward message
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                [Message(DeltaDistribution([1.0, 2.0])), Message(DeltaDistribution([3.0, 4.0])), nothing],
                                DeltaDistribution([4.0, 6.0]))
        # Backward message
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                [nothing, Message(DeltaDistribution([1.0, 2.0])), Message(DeltaDistribution([3.0, 4.0]))],
                                DeltaDistribution([2.0, 2.0]))
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                [Message(DeltaDistribution([1.0, 2.0])), nothing, Message(DeltaDistribution([3.0, 4.0]))],
                                DeltaDistribution([2.0, 2.0]))
    end

    # Tests on Gaussian messages use the update rules from Korl (2005),
    # "A Factor Graph Approach to Signal Modelling, System Identification and Filtering.", Table 4.1.
    context("AdditionNode should propagate a univariate GaussianDistribution") do
        context("AdditionNode should propagate a univariate GaussianDistribution with (m,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0)), nothing],
                                    GaussianDistribution(m=4.0, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=1.0, V=2.0)), Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    [Message(GaussianDistribution(m=1.0, V=2.0)), nothing, Message(GaussianDistribution(m=3.0, V=4.0))],
                                    GaussianDistribution(m=2.0, V=6.0))        
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with (m,W) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=1.0, W=2.0)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    [Message(GaussianDistribution(m=1.0, W=2.0)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))  
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with (xi,V) parametrization") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0)), nothing],
                                    GaussianDistribution(xi=14/6, V=6.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(xi=1.0, V=2.0)), Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    [Message(GaussianDistribution(xi=1.0, V=2.0)), nothing, Message(GaussianDistribution(xi=3.0, V=4.0))],
                                    GaussianDistribution(xi=10/6, V=6.0))  
        end

        context("AdditionNode should propagate a univariate GaussianDistribution with different parametrizations") do
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0)), nothing],
                                    GaussianDistribution(m=4.0, W=4.0/3.0))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=1.0, V=0.5)), Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    [Message(GaussianDistribution(m=1.0, V=0.5)), nothing, Message(GaussianDistribution(m=3.0, W=4.0))],
                                    GaussianDistribution(m=2.0, W=4.0/3.0))  
        end
    end

    context("AdditionNode should propagate a combination of GaussianDistribution and DeltaDistribution") do
        # Forward message
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                [Message(GaussianDistribution(m=1.0, V=0.5)), Message(DeltaDistribution(3.0)), nothing],
                                GaussianDistribution(m=4.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(), 
                                3, 
                                [Message(DeltaDistribution(3.0)), Message(GaussianDistribution(m=1.0, V=0.5)), nothing],
                                GaussianDistribution(m=4.0, V=0.5+tiny))
        #Backward message towards in1
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                [nothing, Message(DeltaDistribution(3.0)), Message(GaussianDistribution(m=1.0, V=0.5))],
                                GaussianDistribution(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(), 
                                1, 
                                [nothing, Message(GaussianDistribution(m=1.0, V=0.5)), Message(DeltaDistribution(3.0))],
                                GaussianDistribution(m=2.0, V=0.5+tiny))
        # Backward message towards in2
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                [Message(DeltaDistribution(3.0)), nothing, Message(GaussianDistribution(m=1.0, V=0.5))],
                                GaussianDistribution(m=-2.0, V=0.5+tiny))
        validateOutboundMessage(AdditionNode(), 
                                2, 
                                [Message(GaussianDistribution(m=1.0, V=0.5)), nothing, Message(DeltaDistribution(3.0))],
                                GaussianDistribution(m=2.0, V=0.5+tiny))
    end

    context("AdditionNode should propagate a multivariate GaussianDistribution") do
        context("AdditionNode should propagate a multivariate GaussianDistribution with (m,V) parametrization") do
            mean = [1.0:3.0]
            variance = reshape([4.0, 3.0, 2.0,
                                3.0, 4.0, 3.0,
                                2.0, 3.0, 4.0], 3, 3)
            # Forward message
            validateOutboundMessage(AdditionNode(), 
                                    3, 
                                    [Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance)), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=mean, V=variance)), Message(GaussianDistribution(m=mean, V=variance))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
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
                                    [Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision)), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(m=mean, W=precision))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
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
                                    [Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance)), nothing],
                                    GaussianDistribution(xi=[1.0, 2.0, 3.0], V=2.0*variance))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(xi=xi, V=variance)), Message(GaussianDistribution(xi=xi, V=variance))],
                                    GaussianDistribution(xi=[0.0, 0.0, 0.0], V=2.0*variance))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
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
                                    [Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision))), nothing],
                                    GaussianDistribution(m=[2.0, 4.0, 6.0], W=0.5*precision))
            # Backward message
            validateOutboundMessage(AdditionNode(), 
                                    1, 
                                    [nothing, Message(GaussianDistribution(m=mean, W=precision)), Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
            validateOutboundMessage(AdditionNode(), 
                                    2, 
                                    [Message(GaussianDistribution(m=mean, W=precision)), nothing, Message(GaussianDistribution(xi=precision*mean, V=inv(precision)))],
                                    GaussianDistribution(m=[0.0, 0.0, 0.0], W=0.5*precision))
        end
    end
end