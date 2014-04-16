module ForneyLab

export Message, Node, Interface, Edge
export calculatemessage!, calculatemessages!, calculateforwardmessage!, calculatebackwardmessage!

abstract Message
abstract Node

type Interface
    # An Interface belongs to a node and is used to send/receive messages.
    # An Interface has exactly one partner interface, with wich it forms an edge.
    # An Interface can be seen as a half-edge, that connects to a partner Interface to form a complete edge.
    # A message from node a to node b is stored at the Interface of node a that connects to an Interface of node b.
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

type Edge
    # An Edge joins two interfaces and has a direction (from tail to head).
    # Edges are mostly useful for code readability, they are not used internally.
    # Forward messages flow in the direction of the Edge (tail to head).
    tail::Interface
    head::Interface

    function Edge(tail::Interface, head::Interface)
        if  typeof(head.message) == Nothing ||
            typeof(tail.message) == Nothing ||
            typeof(head.message) == typeof(tail.message)
            if !is(head.node, tail.node)
                tail.partner = head
                head.partner = tail
                new(tail, head)
            else
                error("Cannot connect two interfaces of the same node: ", typeof(head.node), " ", head.name)
            end
        else
            error("Head and tail message types do not match: ", typeof(head.message), " and ", typeof(tail.message))
        end
    end
end
function Edge(tailNode::Node, headNode::Node)
    # Create an Edge from tailNode to headNode.
    # Use the first free interface on each node.
    tail = nothing
    head = nothing
    for interface in tailNode.interfaces
        if interface.partner==nothing
            tail = interface
            break
        end
    end
    if tail==nothing
        error("Cannot create edge: no free interface on tail node: ", typeof(tailNode), " ", tailNode.name)
    end
    for interface in headNode.interfaces
        if interface.partner==nothing
            head = interface
            break
        end
    end
    if head==nothing
        error("Cannot create edge: no free interface on head node: ", typeof(headNode), " ", headNode.name)
    end

    return Edge(tail, head)
end

# Messages
include("messages.jl")

# Nodes
include("nodes/constant.jl")
include("nodes/multiplication.jl")

#############################
# Generic methods
#############################

function calculatemessage!(interface::Interface, node::Node, messageType::DataType=Message)
    # Calculate the outbound message on a specific interface of a specified node.
    # The message is stored in the specified interface.
    # Optionally, messageType defines the desired type of the calculated message.

    # Sanity check
    if !is(interface.node, node)
        error("Specified interface does not belong to the specified node (", typeof(node), " ", node.name,")")
    end

    # Calculate all inbound messages
    inbound_message_types = Union() # Union of all inbound message types
    for node_interface_id = 1:length(node.interfaces)
        node_interface = node.interfaces[node_interface_id]
        if is(node_interface, interface) continue end
        if node_interface.partner == nothing
            error("Cannot receive messages on disconnected interface ", node_interface_id, " of ", typeof(node), " ", node.name)
        end
        if node_interface.partner.message == nothing
            # Recursive call to calculate required inbound message
            calculatemessage!(node_interface.partner)
            if node_interface.partner.message == nothing
                error("Could not calculate required inbound message on interface ", node_interface_id, " of ", typeof(node), " ", node.name)
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
    calculatemessage!(interface_id, node, inbound_messages, messageType)
end
calculatemessage!(interface::Interface, messageType::DataType=Message) = calculatemessage!(interface, interface.node, messageType)

function calculatemessages!(node::Node)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculatemessage!(interface, node)
    end
end

# Calculate forward/backward messages on an Edge
calculateforwardmessage!(edge::Edge) = calculatemessage!(edge.tail)
calculatebackwardmessage!(edge::Edge) = calculatemessage!(edge.head)

function clearmessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        interface.message = nothing
    end
end

function clearmessages!(edge::Edge)
    # Clear all messages on an edge.
    edge.head.message = nothing
    edge.tail.message = nothing
end

end # module ForneyLab