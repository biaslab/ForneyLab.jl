############################################
# MulitiplicationNode
############################################
# Description:
#   Multiplication with a multiplier that is received over an edge.
#   sink = multiplier * source
#
# Interface ids, (names) and supported message types:
#   1. (source):
#       GaussianMessage
#   2. (multiplier):
#       GeneralMessage{FloatingPoint}
#   3. (sink):
#       GaussianMessage
############################################

export MultiplicationNode

type MultiplicationNode <: Node
    interfaces::Array{Interface,1}
    name::ASCIIString
    multiplier::Interface
    source::Interface
    sink::Interface

    function MultiplicationNode(name::ASCIIString="#undef")
        self = new(Array(Interface, 3), name)
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.source     = self.interfaces[1]
        self.multiplier = self.interfaces[2]
        self.sink       = self.interfaces[3]
        return self
    end
end

function calculatemessage!{T<:Union(GaussianMessage,GeneralMessage)}(
                            interfaceId::Int,
                            node::MultiplicationNode,
                            inboundMessages::Array{T,1},
                            messageType::DataType)
    # TODO: implement the message equations
    node.interfaces[interfaceId].message = GaussianMessage()
end