export ScheduleEntry, Schedule, ExternalSchedule

type ScheduleEntry
    interface::Interface
    message_calculation_rule::Function  # Is called to calculate the message. Default is sumProduct!.
    post_processing::Function           # Optional, a function that performs post-processing on the message. Leave undefined to skip.
    function ScheduleEntry(interface::Interface, message_calculation_rule::Function, post_processing::Union(Nothing,Function)=nothing)
        if post_processing != nothing
            return new(interface, message_calculation_rule, post_processing)
        else
            return new(interface, message_calculation_rule)
        end
    end
end
ScheduleEntry(interface::Interface) = ScheduleEntry(interface, sumProduct!)

typealias Schedule Array{ScheduleEntry, 1}

function convert(::Type{Schedule}, interfaces::Array{Interface, 1}, message_calculation_rule::Function=sumProduct!, post_processing::Union(Nothing,Function)=nothing)
    # Convert a list of interfaces to an actual schedule
    return ScheduleEntry[ScheduleEntry(iface, message_calculation_rule, post_processing) for iface in interfaces]
end

function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule")
    println(io, " #  | {message calc. rule}  | {post-processing}     | {node type} {node name}:{outbound iface index} ({outbound iface name})")
    println(io, " ---------------------------------------------------------------------------------------------------------------------------")
    entry_counter = 1
    for schedule_entry in schedule
        interface = schedule_entry.interface
        msg_calc_rule = schedule_entry.message_calculation_rule
        postproc = (isdefined(schedule_entry, :post_processing)) ? string(schedule_entry.post_processing) : ""
        interface_name = (getName(interface)!="") ? "($(getName(interface)))" : ""

        println(io, " $(entry_counter)$(" "^(3-length(string(entry_counter))))| $(msg_calc_rule)$(" "^(22-length(string(msg_calc_rule))))| $(postproc)$(" "^(22-length(postproc)))| $(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)")
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
