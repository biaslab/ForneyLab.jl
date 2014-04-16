export ConstantNode

# ConstantNode: has just one interface and always sends out a constant message
type ConstantNode <: Node
    constant::Message
    interfaces::Array{Interface,1}
    interface::Interface

    function ConstantNode(constant::Message)
        self = new(constant, Array(Interface, 1))
        # Create interface
        self.interfaces[1] = Interface(self)
        # Init named interface handle
        self.interface = self.interfaces[1]
        return self
    end
end

function calculatemessage{T<:Message}(
                            interfaceId::Int,
                            node::ConstantNode,
                            inboundMessages::Array{T,1},
                            messageType::DataType)
    node.interfaces[interfaceId].message = node.constant
end