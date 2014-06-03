#####################
# Unit tests
#####################

facts("ConstantNode unit tests") do
    context("ConstantNode() should initialize a ConstantNode with 1 interface") do
        node = ConstantNode()
        @fact typeof(node) => ConstantNode
        @fact length(node.interfaces) => 1
        @fact node.out => node.interfaces[1]
    end

    context("ConstantNode should propagate a GaussianMessage") do
        node = ConstantNode(GaussianMessage(m=[2.0], V=[4.0]))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => GaussianMessage
        @fact node.interfaces[1].message.m => [2.0]
        @fact node.interfaces[1].message.V => reshape([4.0], 1, 1)
    end

    context("ConstantNode should propagate a GeneralMessage") do
        node = ConstantNode(GeneralMessage([1.0, 2.0]))
        @fact node.interfaces[1].message => nothing
        msg = ForneyLab.updateNodeMessage!(1, node)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => GeneralMessage
        @fact node.interfaces[1].message.value => [1.0, 2.0]
    end

    context("setValue!() should set the constant message value") do
        node = ConstantNode()
        setValue!(node, GeneralMessage(42))
        @fact node._value.value => 42
    end

    context("getValue() should get the constant message value") do
        node = ConstantNode(GeneralMessage(42))
        @fact getValue(node).value => 42
    end

    context("value should not be directly accessible") do
        node = ConstantNode(GeneralMessage(42))
        @fact_throws(node.value)
    end
end

#####################
# Integration tests
#####################

facts("ConstantNode integration tests") do
    context("pushMessageInvalidations!() should be called in the background by setValue!(node::ConstantNode, value::Message)") do
        (c1, c2, c3, add, equ) = initializeTreeGraph()
        fwd_msg_y = calculateMessage!(equ.interfaces[3])
        # Check message validity after message passing
        # All forward messages should be valid
        @fact c1.out.message_valid => true
        @fact c2.out.message_valid => true
        @fact c3.out.message_valid => true
        @fact add.out.message_valid => true
        @fact equ.interfaces[3].message_valid => true
        # All backward messages should not be valid
        @fact add.in1.message_valid => false
        @fact add.in2.message_valid => false
        @fact equ.interfaces[1].message_valid => false
        @fact equ.interfaces[2].message_valid => false

        # Now update the value of c2 and check the invalidations
        setValue!(c2, GaussianMessage(m=[-1.], V=[2.]))
        # All messages that depend on c2 should be invalid
        @fact c2.out.message_valid => false
        @fact add.in1.message_valid => false
        @fact add.out.message_valid => false
        @fact equ.interfaces[2].message_valid => false
        @fact equ.interfaces[3].message_valid => false
        # Validity of other messages should not be changed
        @fact c1.out.message_valid => true
        @fact c3.out.message_valid => true
        @fact add.in2.message_valid => false
        @fact equ.interfaces[1].message_valid => false
    end
end