module ForneyLab

#
# Messages
#
abstract Message

# Standard Gaussian message.
type GaussianMessage <: Message
    V::Array{Float64}
    m::Array{Float64}
    # TODO: precision and alternative parametrisations
end
GaussianMessage() = GaussianMessage([1.0], [0.0])

# Message is a scalar or vector of arbitrary type.
# Used for passing parameters (such as multiplication factors) to node functions.
type ScalarParameterMessage{T} <: Message
    value::T
end
ScalarParameterMessage() = ScalarParameterMessage(1.0)

#
# Nodes
#
abstract Node

# An interface is a half-edge that carries a message and connects to a pertner interface.
type Interface
    node::Node
    partner::Union(Interface, Nothing)
    message::Union(Message, Nothing)

    # Sanity check for matching message types
    function Interface(node::Node, partner::Interface, message::Message)
        if typeof(partner) == Nothing || typeof(message) == Nothing # Check if message or partner exist
            new(node, partner, message)
        elseif typeof(message) != typeof(partner.message) # Compare message types
            error("Message type of partner does not match with interface message type")
        else
            new(node, partner, message)
        end
    end
end
Interface(node::Node, message::Message) = Interface(node, nothing, message) 
Interface(node::Node) = Interface(node, nothing, nothing)

# interfaces(node) returns an array of interfaces for iterating
function interfaces(node::Node)
    node_interfaces = Array(Interface,0)
    for field in names(node)
        if typeof(getfield(node, field)) == Interface
            node_interfaces = [node_interfaces, getfield(node, field)]
        end
    end

    return node_interfaces
end

# Constant node
type ConstantNode <: Node
    constant::Float64 # TODO: Is this not a property of the interface?
    interface::Interface
    function ConstantNode(constant)
        self = new(constant)
        self.interface = Interface(self)
        return self
    end
end
ConstantNode() = ConstantNode(1.0)

# Mulitiply source with parameter multiplier: sink = multiplier * source
type MultiplicationNode <: Node
    multiplier::Interface
    source::Interface
    sink::Interface
    function MultiplicationNode()
        self = new()
        self.multiplier = Interface(self)
        self.source = Interface(self)
        self.sink = Interface(self)
        return self
    end
end

# Combine two interfaces to form a complete edge.
# type Edge
#     tail::Interface
#     head::Interface

#     function Edge(tail::Interface, head::Interface)
#         if  typeof(head.message) == Nothing ||
#             typeof(tail.message) == Nothing ||
#             typeof(head.message) == typeof(tail.message)
#             tail.partner = head
#             head.partner = tail
#             new(tail, head)
#         else
#             error("Head and tail message types do not match")
#         end
#     end
# end

#
# calculatemessage calculates the outgoing message on the sink interface (half-edge).
# Every combination of node and message types will get an individual signature for the function call.
#
function calculatemessage(
                            interface::Interface,
                            node::ConstantNode,
                            inboundMessageTypes,
                            messageType::DataType)
    print("Message from constant node")
    interface.message = ScalarParameterMessage(node.constant)
end

function calculatemessage{T<:Union(GaussianMessage,ScalarParameterMessage)}(
                            interface::Interface,
                            node::MultiplicationNode,
                            inboundMessageTypes::T,
                            messageType::DataType)
    # Multiplication of a Gaussian message
    #sink_V = multiplierMessage.value' * sourceMessage.V * multiplierMessage.value
    #sink_m = multiplierMessage.value * sourceMessage.m
    print("Message from multiplier node")
    interface.message = GaussianMessage()
end

# Calculate the outbound message (of type messageType) on a specific interface of a specified node
function calculatemessage(interface::Interface, node::Node, messageType::DataType=Message)
    # Sanity check
    if !is(interface.node, node)
        error("Specified interface does not belong to the specified node")
    end

    # Collect all inbound messages
    inbound_message_types = Union() # Union of all inbound message types
    node_interfaces = interfaces(node) # TODO: check if there is a better way to iterate over the interfaces
    for node_interface in node_interfaces
        if is(node_interface, interface) continue end
        if node_interface.partner == nothing
            error("Cannot receive messages on disconnected interface")
        end
        if node_interface.partner.message == nothing
            # Recursive call to calculate required inbound message
            calculatemessage(node_interface.partner)
            if node_interface.partner.message == nothing
                error("Could not calculate required inbound message")
            end
            inbound_message_types = Union(inbound_message_types, typeof(node_interface.partner.message))
        end
    end

    # Calculate the actual message
    print(typeof(inbound_message_types))
    calculatemessage(interface, node, inbound_message_types, messageType)

    # Clear all inbound messages
    for node_interface in node_interfaces
        if is(node_interface, interface) continue end
        node_interface.partner.message = nothing
    end
end

# Calculate the outbound message on a specified interface
calculatemessage(interface::Interface) = calculatemessage(interface, interface.node)

# Calculate the outbound messages on all interfaces of a specific node
function calculatemessage(node::Node)
    for interface in interfaces(node)
        calculatemessage(interface, node)
    end
end

# Calculate forward/backward messages on edge
#calculateforwardmessage(edge::Edge) = calculatemessage(edge.tail)
#calculatebackwardmessage(edge::Edge) = calculatemessage(edge.head)

end # module ForneyLab