export Interface
export handle

"""
An `Interface` belongs to a `FactorNode` and represents a half-edge.
An `Interface` has at most one partner interface, with wich it forms an edge.
"""
mutable struct Interface
    node::FactorNode
    edge::Union{AbstractEdge, Nothing}
    partner::Union{Interface, Nothing}
end

Interface(node::FactorNode) = Interface(node, nothing, nothing)

function show(io::IO, interface::Interface)
    iface_handle = handle(interface)
    (iface_handle == "") || (iface_handle = "($(iface_handle))")
    println(io, "Interface $(findfirst(isequal(interface), interface.node.interfaces)) $(iface_handle) of $(typeof(interface.node)) $(interface.node.id)")
end

"""Return interface handle name"""
function handle(interface::Interface)
    if isdefined(interface.node, :i)
        for h in keys(interface.node.i)
            if (typeof(h)==Symbol || typeof(h)==Int) && (interface.node.i[h] === interface)
                return string(h)
            end
        end
    end

    return ""
end

Base.isless(i1::Interface, i2::Interface) = isless(name(i1), name(i2))

name(iface::Interface) = string(iface.node.id)*handle(iface)
name(::Nothing) = ""

"""
Determines whether interface must be initialized with a breaker message;
i.e. for EP sites, loopy belief propagation, or situations where outbound
messages depend on inbound messages, such as a Nonlinear update without
a given inverse (RTS smoothing).
"""
function requiresBreaker(interface::Interface)
    partner_interface = ultimatePartner(interface)
    (partner_interface == nothing) && return false # Dangling edge
    
    return requiresBreaker(interface, partner_interface, partner_interface.node) # Dispatch to overloaded methods
end
requiresBreaker(interface::Interface, partner_interface::Interface, partner_node::FactorNode) = false # Default, function is overloaded for separate node types
requiresBreaker(::Nothing) = false # Failsafe

"""
Determine the type and dimensionality of the breaker interface message
"""
function breakerParameters(interface::Interface)
    partner_interface = ultimatePartner(interface)
    
    return breakerParameters(interface, partner_interface, partner_interface.node) # Dispatch to overloaded methods
end

"""
Constructs breaker types dictionary for breaker sites
"""
function breakerTypes(breaker_sites::Vector{Interface})
    breaker_types = Dict{Interface, Type}() # Initialize Interface to Message dictionary
    for site in breaker_sites
        (breaker_type, _) = breakerParameters(site)
        breaker_types[site] = breaker_type
    end

    return breaker_types
end