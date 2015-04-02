export ScheduleEntry, Schedule, setPostProcessing!

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

function setPostProcessing!(schedule::Schedule, interface::Interface, post_processing::Function)
    for entry in schedule
        if entry.interface == interface
            entry.post_processing = post_processing # Edit in place
            return
        end
    end
end

function convert(::Type{Schedule}, interfaces::Array{Interface, 1}, message_calculation_rule::Function=sumProduct!, post_processing::Union(Nothing,Function)=nothing)
    # Convert a list of interfaces to an actual schedule
    return ScheduleEntry[ScheduleEntry(iface, message_calculation_rule, post_processing) for iface in interfaces]
end

function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule")
    println(io, "#  | {message calc. func}  | {post-processing}     | {node type} {node name}:{outbound iface index} ({outbound iface name})       |")
    println(io, "---|-----------------------|-----------------------|------------------------------------------------------------------------------|")
    entry_counter = 1
    for schedule_entry in schedule
        interface = schedule_entry.interface
        msg_calc_func = schedule_entry.message_calculation_rule
        postproc = (isdefined(schedule_entry, :post_processing)) ? string(schedule_entry.post_processing) : ""
        interface_name = (name(interface)!="") ? "($(name(interface)))" : ""
        interface_field = "$(typeof(interface.node)) $(interface.node.name):$(findfirst(interface.node.interfaces, interface)) $(interface_name)"
        println(io, "$(pad(string(entry_counter),3))|$(pad(string(msg_calc_func),23))|$(pad(string(postproc),23))|$(pad(interface_field,78))|")
        entry_counter += 1
    end
end

function show(io::IO, nodes::Array{Node, 1})
     # Show node array (possibly an external schedule)
    println(io, "Nodes:")
    for entry in nodes
        println(io, "Node $(entry.name) of type $(typeof(entry))")
    end
end
