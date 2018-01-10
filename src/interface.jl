export Interface
export handle

"""
An Interface belongs to a FactorNode and represents a half-edge.
An Interface has at most one partner interface, with wich it forms an edge.
"""
mutable struct Interface
    node::FactorNode
    edge::Union{AbstractEdge, Void}
    partner::Union{Interface, Void}
end

Interface(node::FactorNode) = Interface(node, nothing, nothing)

function show(io::IO, interface::Interface)
    iface_handle = handle(interface)
    (iface_handle == "") || (iface_handle = "($(iface_handle))")
    println(io, "Interface $(findfirst(interface.node.interfaces, interface)) $(iface_handle) of $(typeof(interface.node)) $(interface.node.id)")
end

"""Return interface handle name"""
function handle(interface::Interface)
    if isdefined(interface.node, :i)
        for h in keys(interface.node.i)
            if (typeof(h)==Symbol || typeof(h)==Int) && is(interface.node.i[h], interface)
                return string(h)
            end
        end
    end

    return ""
end