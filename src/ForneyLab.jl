module ForneyLab

export  Message, Node, CompositeNode, Interface, Edge
export  calculateMessage!, calculateMessages!, calculateForwardMessage!, calculateBackwardMessage!,
        calculateMarginal,
        getMessage, getForwardMessage, getBackwardMessage, setMessage!, setMessageValid!, setForwardMessage!, setBackwardMessage!, clearMessages!,
        generateSchedule, executeSchedule

# Verbosity
verbose = false
setVerbose(verbose_mode=true) = global verbose=verbose_mode
printVerbose(msg) = if verbose println(msg) end

# Helpers
include("helpers.jl")

# Other includes
import Base.show

# Top-level abstracts
abstract Message

abstract Node
show(io::IO, node::Node) = println(io, typeof(node), " with name ", node.name, ".")
abstract CompositeNode <: Node

type Interface
    # An Interface belongs to a node and is used to send/receive messages.
    # An Interface has exactly one partner interface, with wich it forms an edge.
    # An Interface can be seen as a half-edge, that connects to a partner Interface to form a complete edge.
    # A message from node a to node b is stored at the Interface of node a that connects to an Interface of node b.
    node::Node
    partner::Union(Interface, Nothing) # Partner indicates the interface to which it is connected.
    child::Union(Interface, Nothing) # An interface that belongs to a composite has a child, which is the corresponding (effectively the same) interface one lever deeper in the node hierarchy.
    #self_dependent::Bool # If an interface is self-dependent, the outbound message depends on the inbound message on the same interface
    message::Union(Message, Nothing)
    message_valid::Bool # false indicates that message has already been consumed

    # Sanity check for matching message types
    function Interface(node::Node, partner::Union(Interface, Nothing)=nothing, child::Union(Interface, Nothing)=nothing, message::Union(Message, Nothing)=nothing)
        if typeof(partner) == Nothing || typeof(message) == Nothing # Check if message or partner exist
            new(node, partner, child, message, typeof(message)!=Nothing)
        elseif typeof(message) != typeof(partner.message) # Compare message types
            error("Message type of partner does not match with interface message type")
        else
            new(node, partner, child, message, typeof(message)!=Nothing)
        end
    end
end
Interface(node::Node, message::Message) = Interface(node, nothing, nothing, message)
Interface(node::Node) = Interface(node, nothing, nothing, nothing)
show(io::IO, interface::Interface) = println(io, "Interface of $(typeof(interface.node)) with node name $(interface.node.name) holds ", interface.message_valid ? "VALID" : "INVALID", " message of type $(typeof(interface.message)).")
function setMessage!(interface::Interface, message::Message, track_invalidations::Bool=true)
    interface.message = message
    if track_invalidations pushMessageInvalidations!(interface) end
    interface.message_valid = true
end
getMessage(interface::Interface) = interface.message
setMessageValid!(interface::Interface) = interface.message_valid = (interface.message!=nothing)
function show(io::IO, schedule::Array{Interface, 1})
    # Show schedules in a specific way
    println(io, "Message passing schedule:")
    for interface in schedule
        println(io, " $(typeof(interface.node)) $(interface.node.name)")
    end
end

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
                # Partner head and tail, and merge their families
                tail.partner = head
                head.partner = tail
                # Backreferences for tail's children
                child_interface = tail.child
                while child_interface != nothing
                    child_interface.partner = tail.partner
                    child_interface = child_interface.child
                end
                # Backreferences for head's children
                child_interface = head.child
                while child_interface != nothing
                    child_interface.partner = head.partner
                    child_interface = child_interface.child
                end
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
setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
getForwardMessage(edge::Edge) = edge.tail.message
getBackwardMessage(edge::Edge) = edge.head.message

# Messages
include("messages.jl")

# Nodes
include("nodes/addition.jl")
include("nodes/constant.jl")
include("nodes/equality.jl")
include("nodes/fixed_gain.jl")
# Composite nodes
include("nodes/composite/gain_addition.jl")
include("nodes/composite/gain_equality.jl")
include("nodes/composite/general.jl")

#############################
# Generic methods
#############################

function calculateMessage!(outbound_interface::Interface, track_invalidations::Bool=true)
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...")
    schedule = generateSchedule(outbound_interface)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("Executing above schedule...")
    executeSchedule(schedule, track_invalidations)
    printVerbose("calculateMessage!() done.")
    return outbound_interface.message
end

function updateNodeMessage!(outbound_interface::Interface, track_invalidations::Bool=true)
    # Calculate the outbound message based on the inbound messages and the node update function.
    # The resulting message is stored in the specified interface and returned.

    node = outbound_interface.node

    # Determine types of inbound messages
    inbound_message_types = Union() # Union of all inbound message types
    outbound_interface_id = 0
    for node_interface_id = 1:length(node.interfaces)
        node_interface = node.interfaces[node_interface_id]
        if is(node_interface, outbound_interface)
            outbound_interface_id = node_interface_id
            continue
        end
        @assert(node_interface.partner!=nothing, "Cannot receive messages on disconnected interface $node_interface_id of $(typeof(node)) $(node.name)")
        @assert(node_interface.partner.message!=nothing, "There is no inbound message present on interface $node_interface_id of $(typeof(node)) $(node.name)")
        inbound_message_types = Union(inbound_message_types, typeof(node_interface.partner.message))
    end

    # Evaluate node update function
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id")
    msg = updateNodeMessage!(outbound_interface_id, node, inbound_message_types)
    printVerbose(" >> $(msg)")

    # Invalidate everything that depends on the outbound message
    if track_invalidations pushMessageInvalidations!(outbound_interface) end

    # Validate the outbound message
    outbound_interface.message_valid = (typeof(outbound_interface.message)<:Message)

    return msg
end

function calculateMessages!(node::Node, track_invalidations::Bool=true)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculateMessage!(interface, track_invalidations)
    end
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge, track_invalidations::Bool=true) = calculateMessage!(edge.tail, track_invalidations)
calculateBackwardMessage!(edge::Edge, track_invalidations::Bool=true) = calculateMessage!(edge.head, track_invalidations)

function executeSchedule(schedule::Array{Interface, 1}, track_invalidations::Bool=true)
    # Execute a message passing schedule
    for interface in schedule
        updateNodeMessage!(interface, track_invalidations)
    end
    # Return the last message in the schedule
    return schedule[end].message
end

function calculateMarginal(forward_msg::Message, backward_msg::Message)
    # Calculate the marginal from a forward/backward message pair.
    # We calculate the marginal by using an EqualityNode.
    # The forward and backward messages are inbound messages to the EqualityNode.
    # The outbound message is the marginal.
    @assert(typeof(forward_msg)==typeof(backward_msg), "Cannot create marginal from forward/backward messages of different types.")
    eq_node = EqualityNode(3)
    c_node1 = ConstantNode(forward_msg)
    c_node2 = ConstantNode(backward_msg)
    Edge(c_node1.out, eq_node.interfaces[1])
    Edge(c_node2.out, eq_node.interfaces[2])
    marginal_msg = deepcopy(calculateMessage!(eq_node.interfaces[3]))
    return marginal_msg
end

function calculateMarginal(edge::Edge, track_invalidations::Bool=true)
    # Calculate the marginal message on an edge
    # If the required messages are not available, they will be calculated by calling calculateMessage!().
    @assert(typeof(edge.tail)==Interface, "Edge should be bound to a tail interface.")
    @assert(typeof(edge.head)==Interface, "Edge should be bound to a head interface.")
    if edge.tail.message==nothing
        calculateMessage!(edge.tail, track_invalidations)
    end
    if edge.head.message==nothing
        calculateMessage!(edge.head, track_invalidations)
    end
    return calculateMarginal(edge.tail.message, edge.head.message)
end

function clearMessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        interface.message = nothing
        interface.message_valid = false
    end
end

function clearMessages!(edge::Edge)
    # Clear all messages on an edge.
    edge.head.message = nothing
    edge.head.message_valid = false
    edge.tail.message = nothing
    edge.tail.message_valid = false
end

function generateSchedule(outbound_interface::Interface)
    # Generate a schedule that can be executed to calculate the outbound message on outbound_interface.
    # IMPORTANT: the resulting schedule depends on the current (valid) messages stored in the factor graph.
    # The same graph with different messages being (in)valid can and probably will result in a different schedule.
    # When a lot of iterations of the same message passing schedule are required, it can be very beneficial
    # to generate the schedule just once using this function, and then execute the same schedule over and over.
    # This prevents recursions through the graph in every call to calculateMessage!().
    schedule = generateScheduleByDFS(outbound_interface)
end

function generateSchedule(partial_schedule::Array{Interface, 1})
    # Generate a complete schedule based on partial_schedule.
    # A partial schedule only defines the order of a subset of all required messages.
    # This function will find a valid complete schedule that satisfies the partial schedule.
    # IMPORTANT: the resulting schedule depends on the current (valid) messages stored in the factor graph.
    # The same graph with different messages being (in)valid can and probably will result in a different schedule.
    # When a lot of iterations of the same message passing schedule are required, it can be very beneficial
    # to generate the schedule just once using this function, and then execute the same schedule over and over.
    # This prevents recursions through the graph in every call to calculateMessage!().
    schedule = Array(Interface, 0)
    for interface_order_constraint in partial_schedule
        schedule = generateScheduleByDFS(interface_order_constraint, schedule)
    end

    return schedule
end

function generateScheduleByDFS(outbound_interface::Interface, backtrace::Array{Interface, 1}=Array(Interface, 0), call_list::Array{Interface, 1}=Array(Interface, 0))
    # This is a private function that performs a search through the factor graph to generate a schedule.
    # IMPORTANT: the resulting schedule depends on the current (valid) messages stored in the factor graph.
    # This is a recursive implementation of DFS. The recursive calls are stored in call_list.
    # backtrace will hold the backgrace of the recursion.
    node = outbound_interface.node

    # Apply stopping condition for recursion. When the same interface is called twice, this is indicative of an unbroken loop.
    if outbound_interface in call_list
        # Notify the user to break the loop with an initial message
        error("Loop detected around $(outbound_interface) Consider setting an initial message at this interface.")
    else # Stopping condition not reached
        push!(call_list, outbound_interface) # Increment list
    end

    # Check all inbound messages on the other interfaces of the node
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
        if !(node_interface.partner.message_valid) # When the required inbound message is invalid, calculate it anew
            if !(node_interface.partner in backtrace) # Don't recalculate stuff that's already in the schedule.
                # Recursive call
                printVerbose("Recursive call of generateSchedule! on node $(typeof(node_interface.partner.node)) $(node_interface.partner.node.name)")
                generateScheduleByDFS(node_interface.partner, backtrace, call_list)
            end
        end
    end

    # Update call_list and backtrace
    pop!(call_list)
    push!(backtrace, outbound_interface)
end

function pushMessageInvalidations!(interface::Union(Interface, Nothing), interface_is_inbound::Bool=false, call_list::Array{Interface, 1}=Array(Interface, 0))
    # Invalidate all dependencies of a message. We call two messages dependent when one message (parent message) is used for the calculation of the other (child message).
    # Dependence implies that alteration of the parent message invalidates the child message.

    # If interface_is_inbound==false: Invalidate everything that depends on the SENT (outbound) message message on interface.
    # If interface_is_inbound==true:  Invalidate everything that depends on the RECEIVED (inbound) message message on interface.
    # The message on the argument interface itself is NOT INVALIDATED.

    # A call to this function with only an interface as argument expects the interface to be outbound. This is the behaviour the user sees.
    # Internally, consecutive recursive calls to this function are on inbound interfaces from then on, indicated by interface_is_inbound=true.

    # If called with an outbound interface, translate call to inbound interface
    if interface_is_inbound==false
        pushMessageInvalidations!(interface.partner, true)
    end

    # Stopping condition for recursion
    if interface==nothing || interface in call_list
        return
    else
        push!(call_list, interface)
    end

    # Push invalidations through the other interfaces of the node.
    node = interface.node
    for interface_id = 1:length(node.interfaces)
        if is(interface, node.interfaces[interface_id]) continue end # skip the inbound interface
        node.interfaces[interface_id].message_valid = false
        # This performs a DFS through the graph, invalidating all messages that depend on connected_node.interfaces[interface_id].message
        pushMessageInvalidations!(node.interfaces[interface_id].partner, true, call_list)
    end

    # Push invalidations through the internals of the connected node (if applicable)
    if interface.child!=nothing
        pushMessageInvalidations!(interface.child, true, call_list)
    end

    pop!(call_list)
end

function pushMessageInvalidations!(node::Node)
    # Invalidates all outbound messages of node AND all messages that depend on the node's outbound messages.
    for interface in node.interfaces
        interface.message_valid = false
        pushMessageInvalidations!(interface)
    end
end

try
    # Try to load user-defined extensions
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/src/forneylab_extensions.jl")
end

end # module ForneyLab