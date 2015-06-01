export Interface
export clearMessage!, setMessage!, message, handle

type Interface
    # An Interface belongs to a node and is used to send/receive messages.
    # An Interface has exactly one partner interface, with wich it forms an edge.
    # An Interface can be seen as a half-edge, that connects to a partner Interface to form a complete edge.
    # A message from node a to node b is stored at the Interface of node a that connects to an Interface of node b.
    node::Node
    edge::Union(AbstractEdge, Nothing)
    partner::Union(Interface, Nothing)
    message::Union(Message, Nothing)
end
Interface(node::Node) = Interface(node, nothing, nothing, nothing)

function show(io::IO, interface::Interface)
    iface_handle = handle(interface)
    (iface_handle == "") || (iface_handle = "($(iface_handle))")
    println(io, "Interface $(findfirst(interface.node.interfaces, interface)) $(iface_handle) of $(typeof(interface.node)) $(interface.node.id)")
end
function setMessage!(interface::Interface, message::Message)
    interface.message = deepcopy(message)
end
clearMessage!(interface::Interface) = (interface.message=nothing)
message(interface::Interface) = interface.message
function handle(interface::Interface)
    # Return interface handle
    if isdefined(interface.node, :i)
        for h in keys(interface.node.i)
            if typeof(h)==Symbol && is(interface.node.i[h], interface)
                return h
            end
        end
    end
    error("Interface not found")
end

function ensureMessage!{T<:ProbabilityDistribution}(interface::Interface, payload_type::Type{T})
    # Ensure that interface carries a Message{payload_type}, used for in place updates
    if interface.message == nothing || typeof(interface.message.payload) != payload_type
        if payload_type <: DeltaDistribution
            interface.message = Message(payload_type())
        else
            interface.message = Message(vague(payload_type))
        end
    end

    return interface.message
end
