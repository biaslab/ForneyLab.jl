facts("AdditionNode") do
    context("AdditionNode() should initialize an AdditionNode with 3 interfaces") do
        node = AdditionNode()
        @fact typeof(node) => AdditionNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
    end

    context("AdditionNode should add two GeneralMessages") do
        node = AdditionNode()
        inbound_msg_1 = GeneralMessage(2.0)
        inbound_msg_2 = GeneralMessage(3.0)
        # Forward message
        msg = calculateMessage!(3, node, [inbound_msg_1, inbound_msg_2])
        @fact node.interfaces[3].message => msg
        @fact node.interfaces[3].message.value => 5.0
    end

    context("AdditionNode should propagate a univariate GaussianMessage") do
        context("Univariate GaussianMessage with (m,V) parametrization") do
        end

        context("Univariate GaussianMessage with (m,V=0) parametrization") do
        end

        context("Univariate GaussianMessage with (m,W) parametrization") do
        end

        context("Univariate GaussianMessage with (xi,W) parametrization") do
        end
    end

    context("AdditionNode should propagate a multivariate GaussianMessage") do
        context("Multivariate GaussianMessage with (m,V) parametrization") do
        end

        context("Multivariate GaussianMessage with (m,V=0) parametrization") do
        end

        context("Multivariate GaussianMessage with (m,W) parametrization") do
        end

        context("Multivariate GaussianMessage with (xi,W) parametrization") do
        end
    end


end