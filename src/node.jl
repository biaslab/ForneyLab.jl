export Node

abstract Node

Base.deepcopy(::Node) = error("deepcopy(::Node) is not supported. You should use copy(src::Node, id=:new_id) to create an independent copy with no edges attached.")

show(io::IO, node::Node) = println(io, "$(typeof(node)) with id $(node.id)")

function Base.copy(src::Node; id::Symbol = generateNodeId(typeof(src)))
    # Create an independent copy of src with no edges connected to it.
    # The copy is added to the currently active graph.
    # Argument id specifies the id of the copy.

    (!haskey(currentGraph().nodes, id)) || error("The current graph already contains a node with id $(id)")

    # Isolate src from the rest of the graph so deepcopy will only deepcopy src itself (and not other parts of the graph it belongs to)
    src_interface_partners = Dict{Interface,Interface}()
    src_interface_edges = Dict{Interface,Edge}()
    for src_interface in src.interfaces
        src_interface_partners[src_interface] = src_interface.partner
        src_interface.partner = nothing
        src_interface_edges[src_interface] = src_interface.edge
        src_interface.edge = nothing
    end

    # Deepcopy src
    dup = Base.deepcopy_internal(src, ObjectIdDict()) # deepcopy(src) does not work since deepcopy(::Node) is forbidden

    # Reconnect src
    for src_interface in src.interfaces
        src_interface.partner = src_interface_partners[src_interface]
        src_interface.edge = src_interface_edges[src_interface]
    end

    # Update id of copy and add it to the current graph
    dup.id = id
    addNode!(currentGraph(), dup)

    return dup
end

function show(io::IO, nodes::Union{Vector{Node},Set{Node}})
     # Show node array (possibly an external schedule)
    println(io, "Nodes:")
    for node in nodes
        show(io, node)
    end
end
