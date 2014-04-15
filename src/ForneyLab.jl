module ForneyLab

#############################
#   Top-level types
#############################

abstract Message
abstract Node

# Interface: belongs to a node, can be linked to a partner interface and can hold a message.
# It can be viewed as a half-edge that connects to a partner interface to form a complete edge.
type Interface
    node::Node
    partner::Union(Interface, Nothing)
    message::Union(Message, Nothing)
end
Interface(node::Node) = Interface(node, nothing, nothing)

# Edge: combines two interfaces.
# Mostly useful for code readability, not used internally.
# Forward messages are defined to flow from tail to head
type Edge
    tail::Interface
    head::Interface

    function Edge(tail::Interface, head::Interface)
        if  typeof(head.message) == Nothing ||
            typeof(tail.message) == Nothing ||
            typeof(head.message) == typeof(tail.message)
            tail.partner = head
            head.partner = tail
            new(tail, head)
        else
            error("Head and tail message types do not match")
        end
    end
end

#############################
# Messages
#############################

# GaussianMessage: encodes a Gaussian PDF
type GaussianMessage <: Message
    V::Array{Float64}
    m::Array{Float64}
end
GaussianMessage() = GaussianMessage([1.0], [0.0])

# ScalarParameterMessage: simply holds a scalar value
type ScalarParameterMessage{T} <: Message
    value::T
end
ScalarParameterMessage() = ScalarParameterMessage(1.0)

#############################
# Nodes
#############################

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


#############################
# Generic methods
#############################

# Calculate the outbound message (of type messageType) on a specific interface of a specified node
function calculatemessage(interface::Interface, node::Node, messageType::DataType=Message)
    # Sanity check
    if !is(interface.node, node)
        error("Specified interface does not belong to the specified node")
    end

    # Calculate all inbound messages
    inbound_message_types = Union() # Union of all inbound message types
    for node_interface in node.interfaces
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

    # Collect all inbound messages
    inbound_messages = Array(inbound_message_types, length(node.interfaces))
    interface_id = 0
    for node_interface_id = 1:length(node.interfaces)
        node_interface = node.interfaces[node_interface_id]
        if is(node_interface, interface)
            interface_id = node_interface_id
            continue
        end
        inbound_messages[node_interface_id] = node.interfaces[node_interface_id].partner.message
    end

    # Calculate the actual message
    calculatemessage(interface_id, node, inbound_messages, messageType)

    # Clear all inbound messages
    for node_interface in node.interfaces
        if is(node_interface, interface) continue end
        node_interface.partner.message = nothing
    end
end

# Calculate the outbound message on a specified interface
calculatemessage(interface::Interface) = calculatemessage(interface, interface.node)

# Calculate the outbound messages on all interfaces of a specific node
function calculatemessage(node::Node)
    for interface in node.interfaces
        calculatemessage(interface, node)
    end
end

# Calculate forward/backward messages on an edge
calculateforwardmessage(edge::Edge) = calculatemessage(edge.tail)
calculatebackwardmessage(edge::Edge) = calculatemessage(edge.head)

#############################
# Node/message specific methods
#############################
#
# calculatemessage calculates the outgoing message on the sink interface (half-edge).
# Every combination of node and message types will get an individual signature for the function call.
#

function calculatemessage{T<:Message}(
                            interfaceId::Int,
                            node::ConstantNode,
                            inboundMessages::Array{T,1},
                            messageType::DataType)
    node.interfaces[interfaceId].message = node.constant
end

function calculatemessage{T<:Union(GaussianMessage,ScalarParameterMessage)}(
                            interfaceId::Int,
                            node::MultiplicationNode,
                            inboundMessages::Array{T,1},
                            messageType::DataType)
    # Multiplication for Gaussian message
    #sink_V = multiplierMessage.value' * sourceMessage.V * multiplierMessage.value
    #sink_m = multiplierMessage.value * sourceMessage.m
    node.interfaces[interfaceId].message = GaussianMessage()
end

# TODO: move all node/message specific definitions and methods to separate files
end # module ForneyLab