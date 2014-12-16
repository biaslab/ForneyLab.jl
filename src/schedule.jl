export Schedule, ExternalSchedule, SummaryOperation

# Summary operations for calculating outbound messages
# Default is set to SumProduct
typealias SummaryOperation ASCIIString

typealias Schedule Array{(Interface, SummaryOperation), 1}
convert_to_schedule(interfaces::Array{Interface, 1}) = [(intf, "sum_product") for intf in interfaces] # Convert a list of interfaces to an actual schedule

function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule [{node type} {node name}:{outbound iface index} ({outbound iface name})]:")
    for (interface, summary_operation) in schedule
        interface_name = (getName(interface)!="") ? "($(getName(interface)))" : ""
        println(io, " $(summary_operation) update at $(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)")
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