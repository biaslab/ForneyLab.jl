export MultiplicationNode
export calculatemessage

# Mulitiplication node: multiply source with parameter multiplier:
# sink = multiplier * source
type MultiplicationNode <: Node
    interfaces::Array{Interface,1}
    multiplier::Interface
    source::Interface
    sink::Interface

    function MultiplicationNode()
        self = new(Array(Interface, 3))
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

function calculatemessage{T<:Union(GaussianMessage,ScalarParameterMessage)}(
                            interfaceId::Int,
                            node::MultiplicationNode,
                            inboundMessages::Array{T,1},
                            messageType::DataType)
    #sink_V = multiplierMessage.value' * sourceMessage.V * multiplierMessage.value
    #sink_m = multiplierMessage.value * sourceMessage.m
    # TODO: finish this method
    node.interfaces[interfaceId].message = GaussianMessage()
end