############################################
# ConstantNode
############################################
# Description:
#   Simple node with just 1 interface.
#   Always sends out a constant (predefined) message.
#
#           out
#   [value]----->
#
#   out = value
#
# One can pass the value to the constructor, and
# modify the value directly later.
#
# Example:
#   ConstantNode(GaussianDistribution(), name="myconst")
#
# Interface ids, (names) and supported message types:
#   1. (out):
#       Message{Any}
############################################

export ConstantNode

type ConstantNode <: Node
    value::Any
    interfaces::Array{Interface,1}
    name::ASCIIString
    out::Interface

    function ConstantNode(value=1.0; name="unnamed")
        if typeof(value) <: Message
            error("ConstantNode $(name) can not hold value of type Message.")
        end
        self = new(deepcopy(value), Array(Interface, 1), name)
        # Create interface
        self.interfaces[1] = Interface(self)
        # Init named interface handle
        self.out = self.interfaces[1]
        return self
    end
end

# Overload firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::ConstantNode) = (node.out.partner==nothing) ? node.out : error("No free interface on $(typeof(node)) $(node.name)")

function updateNodeMessage!(outbound_interface_id::Int,
                            node::ConstantNode,
                            inbound_messages_value_types::Type{None}=None)
    # Calculate an outbound message. The constant node is the only node that does not accept incoming messages,
    # therefore inbound_messages_value_types is only present for consistency. 
    # This function is not exported, and is only meant for internal use.

    # Just put node.value in the outbound value
    return node.interfaces[outbound_interface_id].message = Message(node.value)
end