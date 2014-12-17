export ScheduleEntry, Schedule, ExternalSchedule

type ScheduleEntry
    interface::Interface
    summary_operation::Symbol # Summary operation for calculating outbound messages. Default is :sumproduct

    function ScheduleEntry(interface::Interface, summary_operation::Symbol)
        summary_operation in [:sumproduct, :sumproduct_sample] || error("Unknown summary operation :$(summary_operation). Please choose between ':sumproduct' and ':sumproduct_sample'.")
        return new(interface, summary_operation)
    end
end
ScheduleEntry(interface::Interface) = ScheduleEntry(interface, :sumproduct)

typealias Schedule Array{ScheduleEntry, 1}
convert(::Type{Schedule}, interfaces::Array{Interface, 1}, summary_operation::Symbol = :sumproduct) = [ScheduleEntry(intf, summary_operation) for intf in interfaces] # Convert a list of interfaces to an actual schedule

function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule")
    println(io, " #  | {summary operation}      | {node type} {node name}:{outbound iface index} ({outbound iface name})")
    println(io, " ------------------------------------------------------------------------------------------------------")
    entry_counter = 1
    for schedule_entry in schedule
        interface = schedule_entry.interface
        summary_operation = schedule_entry.summary_operation
        
        interface_name = (getName(interface)!="") ? "($(getName(interface)))" : ""
        println(io, " $(entry_counter)$(" "^(3-length(string(entry_counter))))| $(summary_operation)$(" "^(25-length(string(summary_operation))))| $(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)")
        entry_counter += 1
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