module ForneyLab

export Message, Node, Interface, Edge
export calculateMessage!, calculateMessages!, calculateForwardMessage!, calculateBackwardMessage!, clearMessages!
export epsilon

#############################
# Helper methods
#############################
verbose = false
epsilon = 1.0e-12 # Resolution for approx. equality checks

# ensureMatrix: ensure that the input is a 2D array or nothing
ensureMatrix{T<:Number}(arr::Array{T, 2}) = arr
ensureMatrix{T<:Number}(arr::Array{T, 1}) = reshape(arr, 1, 1)
ensureMatrix(n::Nothing) = nothing

# isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from.
isRoundedPosDef{T<:FloatingPoint}(arr::Array{T, 2}) = ishermitian(round(arr, 12)) && isposdef(arr, 'L')

# verbosity helpers
setVerbose(verbose_mode=true) = global verbose=verbose_mode
printVerbose(msg) = if verbose println(msg) end

import Base.show

abstract Message

abstract Node
show(io::IO, node::Node) = println(io, typeof(node), " with name ", node.name, ".")

type Interface
    # An Interface belongs to a node and is used to send/receive messages.
    # An Interface has exactly one partner interface, with wich it forms an edge.
    # An Interface can be seen as a half-edge, that connects to a partner Interface to form a complete edge.
    # A message from node a to node b is stored at the Interface of node a that connects to an Interface of node b.
    node::Node
    partner::Union(Interface, Nothing)
    message::Union(Message, Nothing)
    message_valid::Bool # false indicates that message has already been consumed

    # Sanity check for matching message types
    function Interface(node::Node, partner::Union(Interface, Nothing), message::Union(Message, Nothing))
        if typeof(partner) == Nothing || typeof(message) == Nothing # Check if message or partner exist
            new(node, partner, message, typeof(message)!=Nothing)
        elseif typeof(message) != typeof(partner.message) # Compare message types
            error("Message type of partner does not match with interface message type")
        else
            new(node, partner, message, typeof(message)!=Nothing)
        end
    end
end
Interface(node::Node, message::Message) = Interface(node, nothing, message)
Interface(node::Node) = Interface(node, nothing, nothing)
show(io::IO, interface::Interface) = println(io, "Interface of $(typeof(interface.node)) with node name $(interface.node.name) holds ", interface.message_valid ? "VALID" : "INVALID", " message of type $(typeof(interface.message)).")

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
                error("Cannot connect two interfaces of the same node: ", typeof(head.node), " ", head.node.name)
            end
        else
            error("Head and tail message types do not match: ", typeof(head.message), " and ", typeof(tail.message))
        end
    end
end

function Edge(tail_node::Node, head_node::Node)
    # Create an Edge from tailNode to headNode.
    # Use the first free interface on each node.
    tail = nothing
    head = nothing
    for interface in tail_node.interfaces
        if interface.partner==nothing
            tail = interface
            break
        end
    end
    if tail==nothing
        error("Cannot create edge: no free interface on tail node: ", typeof(tail_node), " ", tail_node.name)
    end
    for interface in head_node.interfaces
        if interface.partner==nothing
            head = interface
            break
        end
    end
    if head==nothing
        error("Cannot create edge: no free interface on head node: ", typeof(head_node), " ", head_node.name)
    end

    return Edge(tail, head)
end
show(io::IO, edge::Edge) = println(io, "Edge from ", typeof(edge.tail.node), " with node name ", edge.tail.node.name, " to ", typeof(edge.head.node), " with node name ", edge.head.node.name, ". Forward message type: ", typeof(edge.tail.message), ". Backward message type: ", typeof(edge.head.message), ".")

# Messages
include("messages.jl")

# Nodes
include("nodes/constant.jl")
include("nodes/equality.jl")
include("nodes/fixed_gain.jl")
include("nodes/addition.jl")
# Composite nodes
include("nodes/composite/gain_equality.jl")

#############################
# Generic methods
#############################

function calculateMessage!(outbound_interface::Interface, node::Node)
    # Calculate the outbound message on a specific interface of a specified node.
    # The message is stored in the specified interface.

    # Sanity check
    if !is(outbound_interface.node, node)
        error("Specified interface does not belong to the specified node (", typeof(node), " ", node.name,")")
    end

    # Calculate all inbound messages
    inbound_message_types = Union() # Union of all inbound message types
    outbound_interface_id = 0
    for node_interface_id = 1:length(node.interfaces)
        node_interface = node.interfaces[node_interface_id]
        if is(node_interface, outbound_interface)
            outbound_interface_id = node_interface_id
            continue
        end
        if node_interface.partner == nothing
            error("Cannot receive messages on disconnected interface ", node_interface_id, " of ", typeof(node), " ", node.name)
        end
        if !(node_interface.partner.message_valid)
            # Recursive call to calculate required inbound message
            calculateMessage!(node_interface.partner)
            if !(node_interface.partner.message_valid)
                error("Could not calculate required inbound message on interface ", node_interface_id, " of ", typeof(node), " ", node.name)
            end
            inbound_message_types = Union(inbound_message_types, typeof(node_interface.partner.message))
        end
    end

    # Collect all inbound messages
    inbound_messages = Array(inbound_message_types, length(node.interfaces))
    if (inbound_message_types!=None)
        for node_interface_id = 1:length(node.interfaces)
            if node_interface_id!=outbound_interface_id
                inbound_messages[node_interface_id] = node.interfaces[node_interface_id].partner.message
            end
        end
    end

    # Calculate the actual message
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id")
    msg = updateNodeMessage!(outbound_interface_id, node, inbound_messages)
    printVerbose(" >> $(msg)")

    # Invalidate all consumed inbound messages, validate the outbound message
    for node_interface_id = 1:length(node.interfaces)
        if node_interface_id!=outbound_interface_id
            node.interfaces[node_interface_id].partner.message_valid = false
        end
    end
    outbound_interface.message_valid = (typeof(outbound_interface.message)<:Message)

    return msg
end
calculateMessage!(outbound_interface::Interface) = calculateMessage!(outbound_interface, outbound_interface.node)

function calculateMessages!(node::Node)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculateMessage!(interface, node)
    end
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge) = calculateMessage!(edge.tail)
calculateBackwardMessage!(edge::Edge) = calculateMessage!(edge.head)

function clearMessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        interface.message = nothing
    end
end

function clearMessages!(edge::Edge)
   # Clear all messages on an edge.
   edge.head.message = nothing
   edge.tail.message = nothing
end

end # module ForneyLab