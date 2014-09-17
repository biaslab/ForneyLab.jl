export Schedule, ExternalSchedule

typealias Schedule Array{Interface, 1}
function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule [{node type} {node name}:{outbound iface index} ({outbound iface name})]:")
    for interface in schedule
        interface_name = (getName(interface)!="") ? "($(getName(interface)))" : ""
        println(io, " $(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)")
    end
end

typealias ExternalSchedule Array{Node, 1}
function show(io::IO, nodes::Array{Node, 1})
     # Show node array (possibly an external schedule)
    println(io, "Nodes:")
    for entry in nodes
        println(io, "Node $(entry.name) of type $(typeof(entry))")
    end
end