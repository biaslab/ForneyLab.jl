module ForneyLab

export  Message, Node, CompositeNode, Interface, Schedule, Edge, ExternalSchedule, Message, MessagePayload, ProbabilityDistribution, FactorGraph, Factorization, Subgraph
export  calculateMessage!, calculateMessages!, calculateForwardMessage!, calculateBackwardMessage!,
        calculateMarginal, calculateMarginal!,
        getMessage, getName, getForwardMessage, getBackwardMessage, setMessage!, setMarginal!, setForwardMessage!, setBackwardMessage!, clearMessage!, clearMessages!,
        generateSchedule, executeSchedule, uninformative, getOrCreateMessage, getCurrentGraph, setCurrentGraph, nodes, edges, factorizeMeanField
export  ==
export  current_graph

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
show(io::IO, nodes::Array{Node, 1}) = [show(io, node) for node in nodes]
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
        if isdefined(interface.node, field) && is(getfield(interface.node, field), interface)
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

    function Edge(tail::Interface, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType; add_to_graph::Bool=true)
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

        # Incorporate edge and nodes in current graph
        if add_to_graph # Note: edges connected to clamp nodes are added to the subgraph as well
            graph = getCurrentGraph()
            (length(graph.factorization) == 1) || error("Cannot create Edge in an already factorized graph; first build the graph, then define factorizations.")
            subgraph = graph.factorization[1] # There is only one
            push!(subgraph.internal_edges, self) # Define edge as internal
            push!(subgraph.nodes, tail.node) # Add node to subgraph
            push!(subgraph.nodes, head.node)
        end

        return self
    end
end
Edge(tail::Interface, head::Interface, message_payload_type::DataType; args...) = Edge(tail, head, message_payload_type, message_payload_type; args...)
function Edge(tail::Interface, head::Interface; args...)
    forward_message_payload_type = backward_message_payload_type = GaussianDistribution # Default to Gaussian; we could also do a check whether the nodes accept Gaussians
    (tail.message == nothing) || (forward_message_payload_type = typeof(tail.message.payload))
    (head.message == nothing) || (backward_message_payload_type = typeof(head.message.payload))
    Edge(tail, head, forward_message_payload_type, backward_message_payload_type; args...)
end

# Edge constructors that accept nodes instead of a specific Interface
# firstFreeInterface(node) should be overloaded for nodes with interface-invariant node functions
firstFreeInterface(node::Node) = error("Cannot automatically pick a free interface on non-symmetrical $(typeof(node)) $(node.name)")
Edge(tail_node::Node, head::Interface; args...) = Edge(firstFreeInterface(tail_node), head; args...)
Edge(tail_node::Node, head::Interface, message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), head, message_payload_type; args...)
Edge(tail_node::Node, head::Interface, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), head, forward_message_payload_type, backward_message_payload_type; args...)

Edge(tail::Interface, head_node::Node; args...) = Edge(tail, firstFreeInterface(head_node); args...)
Edge(tail::Interface, head_node::Node, message_payload_type::DataType; args...) = Edge(tail, firstFreeInterface(head_node), message_payload_type; args...)
Edge(tail::Interface, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(tail, firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type; args...)

Edge(tail_node::Node, head_node::Node; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node); args...)
Edge(tail_node::Node, head_node::Node, message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), message_payload_type; args...)
Edge(tail_node::Node, head_node::Node, forward_message_payload_type::DataType, backward_message_payload_type::DataType; args...) = Edge(firstFreeInterface(tail_node), firstFreeInterface(head_node), forward_message_payload_type, backward_message_payload_type; args...)

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

typealias ExternalSchedule Array{Node, 1}
function show(io::IO, schedule::ExternalSchedule)
     # Show external schedule
    println(io, "External schedule:")
    for entry in schedule
        println(io, "Node $(entry.name) of type $(typeof(entry))")
    end
end

setForwardMessage!(edge::Edge, message::Message) = setMessage!(edge.tail, message)
setBackwardMessage!(edge::Edge, message::Message) = setMessage!(edge.head, message)
getForwardMessage(edge::Edge) = edge.tail.message
getBackwardMessage(edge::Edge) = edge.head.message

# FactorGraph and Subgraph types
type Subgraph
    nodes::Set{Node}
    internal_edges::Set{Edge}
    external_edges::Set{Edge}
    internal_schedule::Schedule # Schedule for internal message passing (Dauwels step 2)
    external_schedule::ExternalSchedule # Schedule for updates on nodes connected to external edges (Dauwels step 3)
end
Subgraph() = Subgraph(Set{Node}(), Set{Edge}(), Set{Edge}(), Array(Interface, 0), Array(Node, 0))

type FactorGraph
    factorization::Array{Subgraph, 1} # References to subgraphs
    approximate_marginals::Dict{(Node, Subgraph), MessagePayload} # Approximate margials (q's) at nodes connected to external edges from the perspective of Subgraph
end
# Get and set current graph functions
global current_graph = FactorGraph([Subgraph()], Dict{(Node, Subgraph), MessagePayload}()) # Initialize a current_graph
getCurrentGraph() = current_graph::FactorGraph
setCurrentGraph(graph::FactorGraph) = global current_graph = graph # Set a current_graph

FactorGraph() = setCurrentGraph(FactorGraph([Subgraph()], Dict{(Node, Subgraph), MessagePayload}())) # Initialize a new factor graph; automatically sets current_graph

function conformSubgraph!(subgraph::Subgraph)
    # Updates external edges and nodes field based on internal edges

    subgraph.nodes = Set{Node}()
    subgraph.external_edges = Set{Edge}()

    for internal_edge in subgraph.internal_edges # Collect all nodes in the subgraph
        push!(subgraph.nodes, internal_edge.head.node)
        push!(subgraph.nodes, internal_edge.tail.node)
    end
    subgraph.external_edges = setdiff(edges(subgraph.nodes), subgraph.internal_edges) # External edges are the difference between all edges connected to nodes, and the internal edges

    return subgraph
end

function factorize!(graph::FactorGraph, internal_edges::Set{Edge})
    # Add a subgraph containing the edges specified in internal_edges and conform
    for subgraph in graph.factorization
        setdiff!(subgraph.internal_edges, internal_edges) # Remove edges from existing subgraph
    end
    new_subgraph = Subgraph(Set{Node}(), copy(internal_edges), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
    conformSubgraph!(new_subgraph)

    return new_subgraph
end
factorize!(internal_edges::Set{Edge}) = factorize(getCurrentGraph(), internal_edges)

function factorizeMeanField!(graph::FactorGraph)
    # Generate a mean field factorization
    (length(graph.factorization) == 1) || error("Cannot perform mean field factorization on an already factorized graph.")
    edges = edges(nodes(graph, include_clamps=true, open_composites=false), include_external=false)

    for edge in edges
        if typeof(edge.head.node)==EqualityNode || typeof(edge.tail.node)==EqualityNode
            # Collect all other edges that are connected to this one through equality nodes
            edge_cluster = Set{Edge}()
            connected_edge_buffer = [edge]
            while length(connected_edge_buffer) > 0
                current_edge = pop!(connected_edge_buffer)
                push!(edge_cluster, current_edge)
                for node in [edge.head.node, edge.tail.node]
                    if typoef(node) == EqualityNode
                        for interface in node.interfaces
                            if !is(interface.edge, current_edge) && !(interface.edge in edge_cluster)
                                push!(connected_edge_buffer, interface.edge)
                            end
                        end
                    end
                end
            end
        else
            factorize!(graph, edge)
        end
    end
    return graph
end
factorizeMeanField!() = factorizeMeanField!(getCurrentGraph())

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

function calculateMessage!(outbound_interface::Interface)
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

function pushRequiredInbound!(inbound_array::Array{Any, 1}, outbound_interface::Interface, subgraph::Subgraph)
    # TODO: push required inbound to inbound array

    return inbound_array
end

function updateNodeMessage!(outbound_interface::Interface)
    # Calculate the outbound message based on the inbound messages and the node update function.
    # The resulting message is stored in the specified interface and is returned.
    node = outbound_interface.node

    # inbound_array holds the inbound messages or marginals on every interface of the node (indexed by the interface id)
    inbound_array = Array(Any, 0)
    outbound_interface_id = 0
    for interface_id = 1:length(node.interfaces)
        interface = node.interfaces[interface_id]
        if is(interface, outbound_interface)
            outbound_interface_id = interface_id
        end
        if (!isdefined(outbound_interface, :dependencies) && outbound_interface_id==interface_id) ||
            (isdefined(outbound_interface, :dependencies) && !(interface in outbound_interface.dependencies))
            # Ignore the inbound on this interface
            push!(inbound_array, nothing)
        else # Inbound is required
            inbound_edge = interface.edge
            if interface.partner==nothing
                error("Cannot receive messages on disconnected interface $(interface_id) of $(typeof(node)) $(node.name)")
            elseif interface.partner.message==nothing && inbound_edge.marginal==nothing
                error("There is no inbound message/marginal present on the partner/edge of interface $(interface_id) of $(typeof(node)) $(node.name)")
            end

            # Push the required inbound message or edge/node marginal for the inbound interface to inbound_array
            pushRequiredInbound!(inbound_array, outbound_interface, subgraph)
        end
    end

    # Evaluate node update function
    printVerbose("Calculate outbound message on $(typeof(node)) $(node.name) interface $outbound_interface_id (outbound $(outbound_interface.message_payload_type)):")

    return updateNodeMessage!(node, outbound_interface_id, outbound_interface.message_payload_type, inbound_array...)
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

function executeSchedule(schedule::ExternalSchedule)
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
clearMessages!(graph::FactorGraph) = map(clearMessages!, nodes(graph, open_composites=true))
clearMessages!() = clearMessages!(getCurrentGraph())

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

function addChildNodes!(nodes::Set{Node})
    # Add all child nodes to the nodes set
    composite_nodes_stack = Array(CompositeNode, 0) # Composite nodes to open
    for node in nodes
        if typeof(node) <: CompositeNode
            push!(composite_nodes_stack, node)
        end
    end
    # Open all composite nodes
    while length(composite_nodes_stack) > 0
        composite_node = pop!(composite_nodes_stack)
        for field in names(composite_node)
            if typeof(getfield(composite_node, field)) <: Node
                # Add child
                child_node = getfield(composite_node, field)
                push!(nodes, child_node)
                if typeof(child_node) <: CompositeNode
                    push!(composite_nodes_stack, child_node)
                end
            end
        end
    end
    return nodes
end

function addClampNodes!(nodes::Set{Node})
    # Add all clamp nodes to the nodes set
    for node in nodes
        for interface in node.interfaces
            if typeof(interface.partner.node) == ClampNode
                push!(nodes, interface.partner.node)
            end
        end
    end
    return nodes
end

function nodes(subgraph::Subgraph; open_composites::Bool=true, include_clamps=false)
    all_nodes = copy(subgraph.nodes)

    if open_composites; addChildNodes!(all_nodes); end
    if include_clamps; addClampNodes!(all_nodes); end

    return all_nodes
end

function nodes(graph::FactorGraph; open_composites::Bool=true, include_clamps=false)
    all_nodes = Set{Node}()
    for subgraph in graph.factorization
        union!(all_nodes, subgraph.nodes)
    end

    if open_composites; addChildNodes!(all_nodes); end
    if include_clamps; addClampNodes!(all_nodes); end

    return all_nodes
end
nodes(;args...) = nodes(getCurrentGraph(); args...)

function edges(nodes::Array{Node,1}; include_external=true)
    # Returns the set of edges connected to nodes, including or excluding external edges
    # An external edge has only head or tail in the interfaces belonging to nodes in the nodes array

    edges = Set()
    for node in nodes
        for interface in node.interfaces
            if include_external
                if interface.edge!=nothing && ((interface.edge.tail.node in nodes) || (interface.edge.head.node in nodes))
                    push!(edges, interface.edge)
                end
            else
                if interface.edge!=nothing && (interface.edge.tail.node in nodes) && (interface.edge.head.node in nodes)
                    push!(edges, interface.edge)
                end
            end
        end
    end
    return edges
end

# Utils
include("visualization.jl")

try
    # Try to load user-defined extensions
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/src/forneylab_extensions.jl")
end

end # module ForneyLab