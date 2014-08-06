module ForneyLab

export  Message, Node, CompositeNode, Interface, Schedule, Edge, MarginalSchedule, Message, MessageValue, ProbabilityDistribution
export  calculateMessage!, calculateMessages!, calculateForwardMessage!, calculateBackwardMessage!,
        calculateMarginal, calculateMarginal!,
        getMessage, getName, getForwardMessage, getBackwardMessage, setMessage!, setMarginal!, setForwardMessage!, setBackwardMessage!, clearMessage!, clearMessages!, clearAllMessages!,
        generateSchedule, executeSchedule, uninformative, getOrCreateMessage
export  ==

# Verbosity
verbose = false
setVerbose(verbose_mode=true) = global verbose=verbose_mode
printVerbose(msg) = if verbose println(msg) end

# ForneyLab Helpers
include("helpers.jl")

# Other includes
import Base.show

# Top-level abstracts
abstract AbstractEdge # An Interface belongs to an Edge, but Interface is defined before Edge. Because you can not belong to something undefined, Edge will inherit from AbstractEdge, solving this problem.
abstract ProbabilityDistribution # ProbabilityDistribution can be carried by a Message or an Edge (as marginal)
typealias MessageValue Union(ProbabilityDistribution, Number, Vector, Matrix)
abstract Node
show(io::IO, node::Node) = println(io, typeof(node), " with name ", node.name, ".")
abstract CompositeNode <: Node

# Distributions
include("distributions/gaussian.jl")
include("distributions/gamma.jl")
include("distributions/inverse_gamma.jl")

type Message{T<:MessageValue}
    # Message has a value payload, which must be a subtype of MessageValue.
    value::T
end
Message(value::MessageValue) = Message{typeof(value)}(deepcopy(value))
Message() = Message(1.0)
uninformative(type_in::Type{Float64}) = Message(1.0) # uninformative general message
show(io::IO, message::Message) = println(io, typeof(message), " with value ", message.value, ".")

type Interface
    # An Interface belongs to a node and is used to send/receive messages.
    # An Interface has exactly one partner interface, with wich it forms an edge.
    # An Interface can be seen as a half-edge, that connects to a partner Interface to form a complete edge.
    # A message from node a to node b is stored at the Interface of node a that connects to an Interface of node b.
    node::Node
    edge::Union(AbstractEdge, Nothing)
    partner::Union(Interface, Nothing)  # Partner indicates the interface to which it is connected.
    child::Union(Interface, Nothing)    # An interface that belongs to a composite has a child, which is the corresponding (effectively the same) interface one lever deeper in the node hierarchy.
    message::Union(Message, Nothing)
    message_value_type::DataType        # Indicates the type of the message value that is carried by the interface
                                        # The default is a GaussianDistribution, which can be overwitten by the Edge constructor or setMessage! and setMarginal!
    dependencies::Array{Interface, 1}   # Optional array of interfaces (of the same node) on which the outbound msg on this interface depends.
                                        # If this array is #undef, it means that the outbound msg depends on the inbound msgs on ALL OTHER interfaces of the node.
    internal_schedule::Array{Interface, 1}      # Optional schedule that should be executed to calculate outbound message on this interface.
                                                # The internal_schedule field is used in composite nodes, and holds the schedule for internal message passing.

    # Sanity check for matching message types
    function Interface(node::Node, edge::Union(AbstractEdge, Nothing)=nothing, partner::Union(Interface, Nothing)=nothing, child::Union(Interface, Nothing)=nothing, message::Union(Message, Nothing)=nothing, message_value_type::DataType=GaussianDistribution)
        if typeof(partner) == Nothing || typeof(message) == Nothing # Check if message or partner exist
            return new(node, edge, partner, child, message, message_value_type)
        elseif typeof(message) != typeof(partner.message) # Compare message types
            error("Message type of partner does not match with interface message type")
        else
            return new(node, edge, partner, child, message, message_value_type)
        end
    end
end
Interface(node::Node, message::Message) = Interface(node, nothing, nothing, nothing, message, GaussianDistribution)
Interface(node::Node) = Interface(node, nothing, nothing, nothing, nothing, GaussianDistribution)
show(io::IO, interface::Interface) = println(io, "Interface of $(typeof(interface.node)) with node name $(interface.node.name) holds message of type $(typeof(interface.message)).")
function setMessage!(interface::Interface, message::Message)
    interface.message_value_type = typeof(message.value)
    interface.message = deepcopy(message)
end
clearMessage!(interface::Interface) = (interface.message=nothing)
getMessage(interface::Interface) = interface.message
function getName(interface::Interface)
    # Return interface name
    for field in names(interface.node)
        if is(getfield(interface.node, field), interface)
            return string(field)
        end
    end
    return ""
end

typealias Schedule Array{Interface, 1}
function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule:")
    for interface in schedule
        interface_name = (getName(interface)!="") ? "($(getName(interface)))" : ""
        println(io, " $(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)")
    end
end

type Edge <: AbstractEdge
    # An Edge joins two interfaces and has a direction (from tail to head).
    # Edges are mostly useful for code readability, they are not used internally.
    # Forward messages flow in the direction of the Edge (tail to head).
    # Edges can contain marginals, which are the product of the forward and backward message.

    tail::Interface
    head::Interface
    marginal::Any

    function Edge(tail::Interface, head::Interface, forward_message_value_type::DataType, backward_message_value_type::DataType)
        (forward_message_value_type <: MessageValue) || error("Forward message value type $(forward_message_value_type) not supported")
        (backward_message_value_type <: MessageValue) || error("Backward message value type $(backward_message_value_type) not supported")
        if tail.message != nothing && (typeof(tail.message.value) != forward_message_value_type)
            error("Existing forward message ($(typeof(tail.message))) does not match expected type Message{$(forward_message_value_type)}")
        end
        if head.message != nothing && (typeof(head.message.value) != backward_message_value_type)
            error("Existing backward message ($(typeof(head.message))) does not match expected type Message{$(backward_message_value_type)}")
        end
        (!is(head.node, tail.node)) || error("Cannot connect two interfaces of the same node: ", typeof(head.node), " ", head.node.name)

        self = new(tail, head, nothing)
        # Assign pointed to edge from interfaces
        tail.edge = self
        head.edge = self
        # Partner head and tail, and merge their families
        tail.partner = head
        head.partner = tail
        # Set expected outbound interface message types
        tail.message_value_type = forward_message_value_type
        head.message_value_type = backward_message_value_type

        # Backreferences for tail's children
        child_interface = tail.child
        while child_interface != nothing
            child_interface.partner = tail.partner
            child_interface.edge = self
            child_interface.message_value_type = forward_message_value_type
            child_interface = child_interface.child
        end
        # Backreferences for head's children
        child_interface = head.child
        while child_interface != nothing
            child_interface.partner = head.partner
            child_interface.edge = self
            child_interface.message_value_type = backward_message_value_type
            child_interface = child_interface.child
        end

        return self
    end
end
Edge(tail::Interface, head::Interface, message_value_type::DataType) = Edge(tail, head, message_value_type, message_value_type)
function Edge(tail::Interface, head::Interface)
    forward_message_value_type = backward_message_value_type = GaussianDistribution # Default to Gaussian; we could also do a check whether the nodes accept Gaussians
    (tail.message == nothing) || (forward_message_value_type = typeof(tail.message.value))
    (head.message == nothing) || (backward_message_value_type = typeof(head.message.value))
    Edge(tail, head, forward_message_value_type, backward_message_value_type)
end

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.name)")
Edge(tail_node::Node, head::Interface) = Edge(firstFreeInterface(tail_node), head)
Edge(tail_node::Node, head::Interface, message_value_type::DataType) = Edge(firstFreeInterface(tail_node), head, message_value_type)
Edge(tail_node::Node, head::Interface, forward_message_value_type::DataType, backward_message_value_type::DataType) = Edge(firstFreeInterface(tail_node), head, forward_message_value_type, backward_message_value_type)

Edge(tail::Interface, head_node::Node) = Edge(tail, firstFreeInterface(head_node))
Edge(tail::Interface, head_node::Node, message_value_type::DataType) = Edge(tail, firstFreeInterface(head_node), message_value_type)
Edge(tail::Interface, head_node::Node, forward_message_value_type::DataType, backward_message_value_type::DataType) = Edge(tail, firstFreeInterface(head_node), forward_message_value_type, backward_message_value_type)

Edge(tail_node::Node, head_node::Node) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node))
Edge(tail_node::Node, head_node::Node, message_value_type::DataType) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), message_value_type)
Edge(tail_node::Node, head_node::Node, forward_message_value_type::DataType, backward_message_value_type::DataType) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), forward_message_value_type, backward_message_value_type)

function show(io::IO, edge::Edge)
    println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.name):$(findfirst(edge.tail.node.interfaces, edge.tail)) to $(typeof(edge.head.node)) $(edge.head.node.name):$(findfirst(edge.head.node.interfaces, edge.head)).")
    println(io, "Accepted forward message value type: $(edge.tail.message_value_type).")
    println(io, "Accepted backward message value type: $(edge.head.message_value_type).")
end

function getOrCreateMessage(interface::Interface, assign_value::DataType, arr_dims::Tuple=(1, 1))
    # Looks for a message on interface.
    # When no message is present, it sets and returns a standard message.
    # Otherwise it returns the present message.
    # For Array types we pre-allocate the array size with arr_dims
    if interface.message==nothing
        if assign_value <: ProbabilityDistribution 
            interface.message = Message(assign_value())
        elseif assign_value == Float64
            interface.message = Message(1.0)
        elseif assign_value <: Array{Float64}
            interface.message = Message(zeros(arr_dims))
        else
            error("Unknown assign type argument")
        end
    end
    return interface.message
end
function getOrCreateMarginal(edge::Edge, assign_value::DataType)
    # Looks for a marginal on edge.
    # When no marginal is present, it sets and returns a standard distribution.
    # Otherwise it returns the present marginal. User for fast marginal calculations.
    if edge.marginal==nothing
        if assign_value <: ProbabilityDistribution 
            edge.marginal = assign_value()
        else
            error("Unknown assign type argument")
        end
    end
    return edge.marginal
end

typealias MarginalSchedule Array{Edge, 1}
function show(io::IO, schedule::MarginalSchedule)
    # Show marginal update schedule
    println(io, "Marginal update schedule:")
    for edge in schedule
        println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.name) to $(typeof(edge.head.node)) $(edge.head.node.name)")
    end
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
getForwardMessage(edge::Edge) = edge.tail.message
getBackwardMessage(edge::Edge) = edge.head.message

# Distribution helpers
include("distributions/helpers.jl")
# Nodes
include("nodes/clamp.jl")
include("nodes/addition.jl")
include("nodes/terminal.jl")
include("nodes/equality.jl")
include("nodes/fixed_gain.jl")
include("nodes/gaussian.jl")
# Composite nodes
include("nodes/composite/gain_addition.jl")
include("nodes/composite/gain_equality.jl")
include("nodes/composite/linear.jl")
include("nodes/composite/general.jl")

#############################
# Generic methods
#############################

function calculateMessage!(outbound_interface::Interface)
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...")
    schedule = generateSchedule(outbound_interface)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("Executing above schedule...")
    executeSchedule(schedule)
    printVerbose("calculateMessage!() done.")

    return outbound_interface.message
end

function updateNodeMessage!(outbound_interface::Interface)
    # Calculate the outbound message based on the inbound messages and the node update function.
    # The resulting message is stored in the specified interface and is returned.
    node = outbound_interface.node

    # inbound_array holds the inbound messages or marginals on every interface of the node (indexed by the interface id)
    inbound_array = Array(Union(Message, ProbabilityDistribution, Nothing), length(node.interfaces))
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if is(interface, outbound_interface)
            outbound_interface_id = interface_id
        end
        if (!isdefined(outbound_interface, :dependencies) && outbound_interface_id==interface_id) ||
            (isdefined(outbound_interface, :dependencies) && !(interface in outbound_interface.dependencies))
            # Ignore the inbound on this interface
            inbound_array[interface_id] = nothing
            continue
        end
        if interface.partner==nothing
            error("Cannot receive messages on disconnected interface $(interface_id) of $(typeof(node)) $(node.name)")
        elseif interface.partner.message==nothing
            error("There is no inbound message present on the partner of interface $(interface_id) of $(typeof(node)) $(node.name)")
        end
        # Add message or marginal to the inbound array
        # TODO: PROPERLY CHECK FOR VARIATIONAL
        if interface.edge.marginal!=nothing && (:variational in names(node)) && node.variational
            inbound_array[interface_id] = interface.edge.marginal
        else
            inbound_array[interface_id] = interface.partner.message
        end
    end

    # Evaluate node update function
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id (outbound $(outbound_interface.message_value_type)):")

    return updateNodeMessage!(node, outbound_interface_id, outbound_interface.message_value_type, inbound_array...)
end

function calculateMessages!(node::Node)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculateMessage!(interface)
    end
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge) = calculateMessage!(edge.tail)
calculateBackwardMessage!(edge::Edge) = calculateMessage!(edge.head)

function executeSchedule(schedule::Schedule)
    # Execute a message passing schedule
    for interface in schedule
        updateNodeMessage!(interface)
    end
    # Return the last message in the schedule
    return schedule[end].message
end

function executeSchedule(schedule::MarginalSchedule)
    # Execute a marginal update schedule
    for edge in schedule
        calculateMarginal!(edge)
    end
    # Return the last message in the schedule
    return schedule[end].marginal
end

function setMarginal!(edge::Edge, distribution::Any)
    # Presets the marginal and head- tail messages on edge with the argument message
    # Usually this method is used to set uninformative messages

    # TODO: check this, message value type does not have to be equal to marginal distribution type
    # TODO: marginal should be the product of the message distributions
    edge.tail.message_value_type = edge.head.message_value_type = typeof(distribution)
    edge.head.message = Message(deepcopy(distribution))
    edge.tail.message = Message(deepcopy(distribution))
    edge.marginal = deepcopy(distribution)
end

function calculateMarginal(edge::Edge)
    # Calculates the marginal without writing back to the edge
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
    return calculateMarginal(edge.tail.message.value, edge.head.message.value)
end
function calculateMarginal!(edge::Edge)
    # Calculates and writes the marginal on edge
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
    calculateMarginal!(edge, edge.tail.message.value, edge.head.message.value)
    return edge.marginal
end

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

# Functions to clear ALL MESSAGES in the graph
clearAllMessages!(seed_node::Node) = map(clearMessages!, getAllNodes(seed_node, open_composites=true))
clearAllMessages!(seed_edge::Edge) = map(clearMessages!, getAllNodes(seed_edge.tail, open_composites=true))

function generateSchedule(outbound_interface::Interface)
    # Generate a schedule that can be executed to calculate the outbound message on outbound_interface.
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    # The same graph with different messages being present can (and probably will) result in a different schedule.
    # When a lot of iterations of the same message passing schedule are required, it can be very beneficial
    # to generate the schedule just once using this function, and then execute the same schedule over and over.
    # This prevents having to generate the same schedule in every call to calculateMessage!().
    schedule = generateScheduleByDFS(outbound_interface)
end

function generateSchedule(partial_schedule::Schedule)
    # Generate a complete schedule based on partial_schedule.
    # A partial schedule only defines the order of a subset of all required messages.
    # This function will find a valid complete schedule that satisfies the partial schedule.
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    # The same graph with different messages being present can (and probably will) result in a different schedule.
    # When a lot of iterations of the same message passing schedule are required, it can be very beneficial
    # to generate the schedule just once using this function, and then execute the same schedule over and over.
    # This prevents having to generate the same schedule in every call to calculateMessage!().
    schedule = Array(Interface, 0)
    for interface_order_constraint in partial_schedule
        schedule = generateScheduleByDFS(interface_order_constraint, schedule)
    end

    return schedule
end

function generateScheduleByDFS(outbound_interface::Interface, backtrace::Schedule=Array(Interface, 0), call_list::Array{Interface, 1}=Array(Interface, 0))
    # This is a private function that performs a search through the factor graph to generate a schedule.
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    # This is a recursive implementation of DFS. The recursive calls are stored in call_list.
    # backtrace will hold the backtrace.
    node = outbound_interface.node

    # Apply stopping condition for recursion. When the same interface is called twice, this is indicative of an unbroken loop.
    if outbound_interface in call_list
        # Notify the user to break the loop with an initial message
        error("Loop detected around $(outbound_interface) Consider setting an initial message somewhere in this loop.")
    else # Stopping condition not satisfied
        push!(call_list, outbound_interface)
    end

    # Check all inbound messages on the other interfaces of the node
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if is(interface, outbound_interface)
            outbound_interface_id = interface_id
        end
        if (!isdefined(outbound_interface, :dependencies) && outbound_interface_id==interface_id) ||
           (isdefined(outbound_interface, :dependencies) && !(interface in outbound_interface.dependencies))
            continue
        end
        @assert(interface.partner!=nothing, "Disconnected interface should be connected: interface #$(interface_id) of $(typeof(node)) $(node.name)")

        if interface.partner.message == nothing # Required message missing.
            if !(interface.partner in backtrace) # Don't recalculate stuff that's already in the schedule.
                # Recursive call
                printVerbose("Recursive call of generateSchedule! on node $(typeof(interface.partner.node)) $(interface.partner.node.name)")
                generateScheduleByDFS(interface.partner, backtrace, call_list)
            end
        end
    end

    # Update call_list and backtrace
    pop!(call_list)

    return push!(backtrace, outbound_interface)
end

function getAllNodes(seed_node::Node, node_array::Array{Node,1}=Array(Node,0); open_composites::Bool=true)
    # Return a list of all nodes
    # If open_composites is false, the search will not pass through composite node boundaries.
    if !(seed_node in node_array)
        push!(node_array, seed_node)
        # Partners
        for interface in seed_node.interfaces
            if interface.partner != nothing
                if interface.partner.partner==interface || open_composites
                    # Recursion
                    getAllNodes(interface.partner.node, node_array, open_composites=open_composites)
                end
            end
        end
        # Children
        if open_composites && typeof(seed_node)<:CompositeNode
            for field in names(seed_node)
                if typeof(getfield(seed_node, field)) <: Node
                    # Recursion
                    getAllNodes(getfield(seed_node, field), node_array, open_composites=open_composites)
                end
            end
        end        
    end

    return node_array
end

# Utils
include("visualization.jl")

try
    # Try to load user-defined extensions
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/src/forneylab_extensions.jl")
end

end # module ForneyLab