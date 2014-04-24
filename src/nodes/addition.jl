############################################
# AdditionNode
############################################
# Description:
#   Addition of two messages of same type.
#   out = in1 + in2
#   Example:
#       AdditionNode(; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       GaussianMessage
#       GeneralMessage
#   2. (in2):
#       GaussianMessage
#       GeneralMessage
#   3. (out):
#       GaussianMessage
#       GeneralMessage
############################################

export AdditionNode

type AdditionNode <: Node
    name::ASCIIString
    interfaces::Array{Interface,1}
    in1::Interface
    in2::Interface
    out::Interface

    function AdditionNode(;args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        self = new(name, Array(Interface, 3))
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
        self.out = self.interfaces[3]
        return self
    end
end

# ############################################
# # GeneralMessage methods
# ############################################

# Calculations for a general message type
function calculateMessage!( outbound_interface_id::Int,
                            node::AdditionNode,
                            inbound_messages::Array{GeneralMessage,1})
    if outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward message
        msg = GeneralMessage()
    elseif outbound_interface_id == 3
        # Forward message
        msg = GeneralMessage(inbound_messages[1].value + inbound_messages[2].value)
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg
end