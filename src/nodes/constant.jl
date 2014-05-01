############################################
# ConstantNode
############################################
# Description:
#   Simple node with just 1 interface.
#   Always sends out a constant (predefined) message.
#
# Interface ids, (names) and supported message types:
#   1. (interface):
#       Message
############################################

export ConstantNode

type ConstantNode <: Node
    constant::Message
    interfaces::Array{Interface,1}
    name::ASCIIString
    interface::Interface

    function ConstantNode(constant::Message; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        self = new(constant, Array(Interface, 1), name)
        # Create interface
        self.interfaces[1] = Interface(self)
        # Init named interface handle
        self.interface = self.interfaces[1]
        return self
    end
end
ConstantNode(; args...) = ConstantNode(GeneralMessage(1.0); args...)

function updateNodeMessage!(outbound_interface_id::Int,
                            node::ConstantNode,
                            inbound_messages::Array{None, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # Just pass the unaltered message through
    return node.interfaces[outbound_interface_id].message = deepcopy(node.constant)
end