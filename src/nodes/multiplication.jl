############################################
# MulitiplicationNode
############################################
# Description:
#   Multiplication of two messages.
#   out = in1 * in2
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

export MultiplicationNode

type MultiplicationNode <: Node
    interfaces::Array{Interface,1}
    name::ASCIIString
    in1::Interface
    in2::Interface
    out::Interface

    function MultiplicationNode(;args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        self = new(Array(Interface, 3), name)
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

function calculatemessage!{T<:Union(GaussianMessage,GeneralMessage)}(
                            outbound_interface_id::Int,
                            node::MultiplicationNode,
                            inbound_messages::Array{T,1})
    if outbound_interface_id == 1 # message to source interface

    elseif outbound_interface_id == 2 # message to source interface

    elseif outbound_interface_id == 3 # message to source interface

    end
    node.interfaces[outbound_interface_id].message = GaussianMessage()
end