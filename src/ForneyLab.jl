module ForneyLab

export  Message, Node, CompositeNode, Interface, Schedule, Edge, MarginalSchedule, Message, MessagePayload, ProbabilityDistribution, Factorization
export  calculateMessage!, calculateMessages!, calculateForwardMessage!, calculateBackwardMessage!,
        calculateMarginal, calculateMarginal!,
        getMessage, getName, getForwardMessage, getBackwardMessage, setMessage!, setMarginal!, setForwardMessage!, setBackwardMessage!, clearMessage!, clearMessages!, clearAllMessages!,
        generateSchedule, executeSchedule, uninformative, getOrCreateMessage, getLocalFactorization
export  ==

# Verbosity
verbose = false
setVerbose(verbose_mode=true) = global verbose=verbose_mode
printVerbose(msg) = if verbose println(msg) end

# ForneyLab helpers
include("helpers.jl")

# Other includes
import Base.show

# Top-level abstracts
abstract AbstractEdge # An Interface belongs to an Edge, but Interface is defined before Edge. Because you can not belong to something undefined, Edge will inherit from AbstractEdge, solving this problem.
abstract ProbabilityDistribution # ProbabilityDistribution can be carried by a Message or an Edge (as marginal)
typealias MessagePayload Union(ProbabilityDistribution, Number, Vector, Matrix)
abstract Node
show(io::IO, node::Node) = println(io, typeof(node), " with name ", node.name, ".")
abstract CompositeNode <: Node

# Distributions
include("distributions/gaussian.jl")
include("distributions/gamma.jl")
include("distributions/inverse_gamma.jl")
include("distributions/normal_gamma.jl")
include("distributions/students_t.jl")

type Message{T<:MessagePayload}
    # Message has a payload, which must be a subtype of MessagePayload.
    payload::T
end
Message(payload::MessagePayload) = Message{typeof(payload)}(deepcopy(payload))
Message() = Message(1.0)
==(msg1::Message, msg2::Message) = msg1.payload == msg2.payload
show(io::IO, message::Message) = println(io, typeof(message), " with payload ", message.payload, ".")

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
    message_payload_type::DataType      # Indicates the type of the message payload that is carried by the interface
                                        # The default is a GaussianDistribution, which can be overwitten by the Edge constructor or setMessage! and setMarginal!
    dependencies::Array{Interface, 1}   # Optional array of interfaces (of the same node) on which the outbound msg on this interface depends.
                                        # If this array is #undef, it means that the outbound msg depends on the inbound msgs on ALL OTHER interfaces of the node.
    internal_schedule::Array{Interface, 1}      # Optional schedule that should be executed to calculate outbound message on this interface.
                                                # The internal_schedule field is used in composite nodes, and holds the schedule for internal message passing.

    # Sanity check for matching message types
    function Interface(node::Node, edge::Union(AbstractEdge, Nothing)=nothing, partner::Union(Interface, Nothing)=nothing, child::Union(Interface, Nothing)=nothing, message::Union(Message, Nothing)=nothing, message_payload_type::DataType=GaussianDistribution)
        if typeof(partner) == Nothing || typeof(message) == Nothing # Check if message or partner exist
            return new(node, edge, partner, child, message, message_payload_type)
        elseif typeof(message) != typeof(partner.message) # Compare message types
            error("Message type of partner does not match with interface message type")
        else
            return new(node, edge, partner, child, message, message_payload_type)
        end
    end
end
Interface(node::Node, message::Message) = Interface(node, nothing, nothing, nothing, message, GaussianDistribution)
Interface(node::Node) = Interface(node, nothing, nothing, nothing, nothing, GaussianDistribution)
show(io::IO, interface::Interface) = println(io, "Interface of $(typeof(interface.node)) with node name $(interface.node.name) holds message of type $(typeof(interface.message)).")
function setMessage!(interface::Interface, message::Message)
    interface.message_payload_type = typeof(message.payload)
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

    function Edge(tail::Interface, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType)
        (forward_message_payload_type <: MessagePayload) || error("Forward message payload type $(forward_message_payload_type) not supported")
        (backward_message_payload_type <: MessagePayload) || error("Backward message payload type $(backward_message_payload_type) not supported")
        if tail.message != nothing && (typeof(tail.message.payload) != forward_message_payload_type)
            error("Existing forward message ($(typeof(tail.message))) does not match expected type Message{$(forward_message_payload_type)}")
        end
        if head.message != nothing && (typeof(head.message.payload) != backward_message_payload_type)
            error("Existing backward message ($(typeof(head.message))) does not match expected type Message{$(backward_message_payload_type)}")
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
        tail.message_payload_type = forward_message_payload_type
        head.message_payload_type = backward_message_payload_type

        # Backreferences for tail's children
        child_interface = tail.child
        while child_interface != nothing
            child_interface.partner = tail.partner
            child_interface.edge = self
            child_interface.message_payload_type = forward_message_payload_type
            child_interface = child_interface.child
        end
        # Backreferences for head's children
        child_interface = head.child
        while child_interface != nothing
            child_interface.partner = head.partner
            child_interface.edge = self
            child_interface.message_payload_type = backward_message_payload_type
            child_interface = child_interface.child
        end

        return self
    end
end
Edge(tail::Interface, head::Interface, message_payload_type::DataType) = Edge(tail, head, message_payload_type, message_payload_type)
function Edge(tail::Interface, head::Interface)
    forward_message_payload_type = backward_message_payload_type = GaussianDistribution # Default to Gaussian; we could also do a check whether the nodes accept Gaussians
    (tail.message == nothing) || (forward_message_payload_type = typeof(tail.message.payload))
    (head.message == nothing) || (backward_message_payload_type = typeof(head.message.payload))
    Edge(tail, head, forward_message_payload_type, backward_message_payload_type)
end

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.name)")
Edge(tail_node::Node, head::Interface) = Edge(firstFreeInterface(tail_node), head)
Edge(tail_node::Node, head::Interface, message_payload_type::DataType) = Edge(firstFreeInterface(tail_node), head, message_payload_type)
Edge(tail_node::Node, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType) = Edge(firstFreeInterface(tail_node), head, forward_message_payload_type, backward_message_payload_type)

Edge(tail::Interface, head_node::Node) = Edge(tail, firstFreeInterface(head_node))
Edge(tail::Interface, head_node::Node, message_payload_type::DataType) = Edge(tail, firstFreeInterface(head_node), message_payload_type)
Edge(tail::Interface, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType) = Edge(tail, firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type)

Edge(tail_node::Node, head_node::Node) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node))
Edge(tail_node::Node, head_node::Node, message_payload_type::DataType) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), message_payload_type)
Edge(tail_node::Node, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type)

function show(io::IO, edge::Edge)
    println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.name):$(findfirst(edge.tail.node.interfaces, edge.tail)) to $(typeof(edge.head.node)) $(edge.head.node.name):$(findfirst(edge.head.node.interfaces, edge.head)).")
    println(io, "Accepted forward message payload type: $(edge.tail.message_payload_type).")
    println(io, "Accepted backward message payload type: $(edge.head.message_payload_type).")
end

function getOrCreateMessage(interface::Interface, assign_payload::DataType, arr_dims::Tuple=(1, 1))
    # Looks for a message on interface.
    # When no message is present, it sets and returns a standard message.
    # Otherwise it returns the present message.
    # For Array types we pre-allocate the array size with arr_dims
    if interface.message==nothing
        if assign_payload <: ProbabilityDistribution 
            interface.message = Message(assign_payload())
        elseif assign_payload == Float64
            interface.message = Message(1.0)
        elseif assign_payload <: Array{Float64}
            interface.message = Message(zeros(arr_dims))
        else
            error("Unknown assign type argument")
        end
    end
    return interface.message
end
function getOrCreateMarginal(edge::Edge, assign_payload::DataType)
    # Looks for a marginal on edge.
    # When no marginal is present, it sets and returns a standard distribution.
    # Otherwise it returns the present marginal. Used for fast marginal calculations.
    if edge.marginal==nothing
        if assign_payload <: ProbabilityDistribution 
            edge.marginal = assign_payload()
        else
            error("Unknown assign type argument")
        end
    end
    return edge.marginal
end
function getOrCreateMarginal(node::Node, assign_payload::DataType)
    # Looks for a marginal on node.
    # When no marginal is present, it sets and returns a standard distribution.
    # Otherwise it returns the present marginal. Used for fast marginal calculations.
    if node.marginal==nothing
        if assign_payload <: ProbabilityDistribution 
            node.marginal = assign_payload()
        else
            error("Unknown assign type argument")
        end
    end
    return node.marginal
end

typealias MarginalSchedule{T<:Union(Edge, Node)} Array{T, 1}
function show(io::IO, schedule::MarginalSchedule)
    # Show marginal update schedule
    println(io, "Marginal update schedule:")
    for entry in schedule
        if typeof(entry) == Edge
            println(io, "Edge from $(typeof(edge.tail.node)) $(edge.tail.node.name) to $(typeof(edge.head.node)) $(edge.head.node.name)")
        else # Node type
            println(io, "Node $(entry.name) of type $(typeof(entry))")
        end
    end
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
getForwardMessage(edge::Edge) = edge.tail.message
getBackwardMessage(edge::Edge) = edge.head.message

# A Factorization is a dictionary that defines the full factorization for the graph of the approximating (q) distribution.
# It couples an integer value to an edge, indicating that this edge is an internal edge of the specified subgraph (Dauwels, 2007).
typealias Factorization Dict{Edge, Int}
function getLocalFactorization(node::Node, factorization::Factorization)
    # Returns a list of the subgraph ids for the edges around 'node' that are defined in 'factorization'
    list = Array(Int64, 0)
    for interface = node.interfaces
        if haskey(factorization, interface.edge)
            push!(list, factorization[interface.edge])
        end
    end
    return list
end

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
# Methods for calculating marginals
include("distributions/calculate_marginal.jl")


#############################
# Generic methods
#############################

function calculateMessage!(outbound_interface::Interface, factorization::Factorization)
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...")
    schedule = generateSchedule(outbound_interface)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("Executing above schedule...")
    executeSchedule(schedule, factorization)
    printVerbose("calculateMessage!() done.")

    return outbound_interface.message
end
calculateMessage!(outbound_interface::Interface) = calculateMessage!(outbound_interface, Factorization())

function insertRequiredInbound!(inbound_array::Array{Union(Message, MessagePayload, Nothing), 1},
                       node::Node,
                       factorization::Factorization,
                       inbound_edge::Edge,
                       outbound_interface::Interface,
                       joint_set_subgraph_id::Array{Int, 1},
                       inbound_interface::Interface,
                       inbound_interface_id::Int)
    # Used by updateNodeMessage(Interface, Factorization)
    # Uses the graph factorization and local node context to determine the required inbound type.
    # Depending on the context, the required inbound can be:
    # - a message from the inbound partner interface,
    # - a marginal from the inbound edge,
    # - a marginal from the node.

    # Collect the incoming edges and marginals
    if !(isdefined(inbound_array, inbound_interface_id) && inbound_array[inbound_interface_id] == nothing) # (Double) check if inbound is required
        if !haskey(factorization, inbound_edge)
            # If inbound_edge is no key in factorization, then take the message from the inbound interface (standard SP)
            inbound_array[inbound_interface_id] = inbound_interface.partner.message
        else # Edge is marked as internal edge on a subgraph
            if factorization[inbound_edge] == factorization[outbound_interface.edge]
                # If the outbound and the inbound edge are in the same subgraph, take the message from the inbound's partner interface
                inbound_array[inbound_interface_id] = inbound_interface.partner.message
            else # The inbound and outbound edges are not in the same subgraph, we need to take a node or edge marginal
                if factorization[inbound_edge] in joint_set_subgraph_id
                    # The inbound edge is part of a joint factorization, take the node marginal.
                    inbound_array[inbound_interface_id] = node.marginal
                else # The inbound edge is part of a univariate factorization, take the inbound edge marginal
                    inbound_array[inbound_interface_id] = inbound_edge.marginal
                end
            end
        end
    end

    return inbound_array
end

function updateNodeMessage!(outbound_interface::Interface, factorization::Factorization)
    # Calculate the outbound message based on the inbound messages and the node update function.
    # The resulting message is stored in the specified interface and is returned.
    node = outbound_interface.node

    local_edge_subgraph_ids = getLocalFactorization(node, factorization) # Subgraph id's for edges around 'node'
    joint_set_subgraph_id = getDuplicatedIds(local_edge_subgraph_ids) # Find the id's that are present more than once (these are the joint factorizations)
    (length(joint_set_subgraph_id) <= 1) || error("There is more than one joint factorization on $(typeof(node)) $(node.name). This is not supported at the moment because a node can only hold one marginal.")

    # inbound_array holds the inbound messages or marginals on every interface of the node (indexed by the interface id)
    inbound_array = Array(Union(Message, MessagePayload, Nothing), length(node.interfaces))
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
        inbound_edge = interface.edge
        if interface.partner==nothing
            error("Cannot receive messages on disconnected interface $(interface_id) of $(typeof(node)) $(node.name)")
        elseif interface.partner.message==nothing && inbound_edge.marginal==nothing
            error("There is no inbound message/marginal present on the partner/edge of interface $(interface_id) of $(typeof(node)) $(node.name)")
        end

        # Plug in the required inbound message or edge/node marginal for the inbound interface in inbound_array
        insertRequiredInbound!(inbound_array, node, factorization, inbound_edge, outbound_interface, joint_set_subgraph_id, interface, interface_id)
    end

    # Evaluate node update function
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id (outbound $(outbound_interface.message_payload_type)):")

    return updateNodeMessage!(node, outbound_interface_id, outbound_interface.message_payload_type, inbound_array...)
end
updateNodeMessage!(outbound_interface::Interface) = updateNodeMessage!(outbound_interface, Factorization())

function calculateMessages!(node::Node, factorization::Factorization)
    # Calculate the outbound messages on all interfaces of node.
    for interface in node.interfaces
        calculateMessage!(interface, factorization)
    end
end
calculateMessages!(node::Node) = calculateMessages!(node, Factorization())

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge, factorization::Factorization) = calculateMessage!(edge.tail, factorization)
calculateForwardMessage!(edge::Edge) = calculateForwardMessage!(edge, Factorization())
calculateBackwardMessage!(edge::Edge, factorization::Factorization) = calculateMessage!(edge.head, factorization)
calculateBackwardMessage!(edge::Edge) = calculateBackwardMessage!(edge, Factorization())

function executeSchedule(schedule::Schedule, factorization::Factorization)
    # Execute a message passing schedule
    for interface in schedule
        updateNodeMessage!(interface, factorization)
    end
    # Return the last message in the schedule
    return schedule[end].message
end
executeSchedule(schedule::Schedule) = executeSchedule(schedule, Factorization())

function executeSchedule(schedule::MarginalSchedule)
    # Execute a marginal update schedule
    for entry in schedule
        calculateMarginal!(entry)
    end
    # Return the last message in the schedule
    return schedule[end].marginal
end

function setMarginal!(edge::Edge, distribution::Any)
    # Presets the marginal on edge with the argument message
    # Usually this method is used to set an uninformative distribution
    edge.marginal = deepcopy(distribution)
end
function setMarginal!(node::Node, distribution::Any)
    # Presets the marginal on a node
    # Usually this method is used to set an uninformative distribution
    node.marginal = deepcopy(distribution)
end

function calculateMarginal(edge::Edge)
    # Calculates the marginal without writing back to the edge
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
    return calculateMarginal(edge.tail.message.payload, edge.head.message.payload)
end
function calculateMarginal!(edge::Edge)
    # Calculates and writes the marginal on edge
    @assert(edge.tail.message != nothing, "Edge should hold a forward message.")
    @assert(edge.head.message != nothing, "Edge should hold a backward message.")
    calculateMarginal!(edge, edge.tail.message.payload, edge.head.message.payload)
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

function getAllNodes(seed_node::Node, node_array::Array{Node,1}=Array(Node,0); open_composites::Bool=true, include_clamps=false)
    # Return a list of all nodes
    # If open_composites is false, the search will not pass through composite node boundaries.
    if include_clamps==false && typeof(seed_node)==ClampNode
        return node_array
    end
    if !(seed_node in node_array)
        push!(node_array, seed_node)
        # Partners
        for interface in seed_node.interfaces
            if interface.partner != nothing
                if interface.partner.partner==interface || open_composites
                    # Recursion
                    getAllNodes(interface.partner.node, node_array, open_composites=open_composites, include_clamps=include_clamps)
                end
            end
        end
        # Children
        if open_composites && typeof(seed_node)<:CompositeNode
            for field in names(seed_node)
                if typeof(getfield(seed_node, field)) <: Node
                    # Recursion
                    getAllNodes(getfield(seed_node, field), node_array, open_composites=open_composites, include_clamps=include_clamps)
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