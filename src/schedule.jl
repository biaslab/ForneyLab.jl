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

Base.deepcopy(::ScheduleEntry) = error("deepcopy(::ScheduleEntry) is not possible. You should construct a new ScheduleEntry or use copy(::ScheduleEntry).")

Base.copy(src::ScheduleEntry) = ScheduleEntry(src.interface, src.message_calculation_rule, isdefined(src.post_processing) ? src.post_processing : nothing)

function setPostProcessing!(schedule_entry::ScheduleEntry, post_processing::Function)
    schedule_entry.post_processing = post_processing
    return schedule_entry
end

typealias Schedule Array{ScheduleEntry, 1}

Base.deepcopy(src::Schedule) = ScheduleEntry[copy(entry) for entry in src]

function setPostProcessing!(schedule::Schedule, interface::Interface, post_processing::Function)
    for entry in schedule
        if entry.interface == interface
            entry.post_processing = post_processing # Edit in place
            return
        end
    end
end

# Convert interfaces to schedule
convert(::Type{ScheduleEntry}, interface::Interface) = ScheduleEntry(interface, sumProduct!) # Default assumes conversion to sum product
convert(::Type{ScheduleEntry}, interface::Interface, message_calculation_rule::Function) = ScheduleEntry(interface, message_calculation_rule)
convert(::Type{Schedule}, interfaces::Array{Interface, 1}, message_calculation_rule::Function) = ScheduleEntry[convert(ScheduleEntry, iface, message_calculation_rule) for iface in interfaces]
convert(::Type{ScheduleEntry}, interface::Interface, message_calculation_rule::Function, post_proc::Function) = ScheduleEntry(interface, message_calculation_rule, post_proc)
convert(::Type{Schedule}, interfaces::Array{Interface, 1}, message_calculation_rule::Function, post_proc::Function) = ScheduleEntry[convert(ScheduleEntry, iface, message_calculation_rule, post_proc) for iface in interfaces]

function show(io::IO, schedule::Schedule)
    # Show schedules in a specific way
    println(io, "Message passing schedule (entry: node [interface], rule)")
    println(io, "------------------------------------------------------")
    entry_counter = 1
    for schedule_entry in schedule
        interface = schedule_entry.interface
        msg_calc_func = schedule_entry.message_calculation_rule
        postproc = (isdefined(schedule_entry, :post_processing)) ? string(schedule_entry.post_processing) : ""
        interface_handle = (handle(interface)!="") ? "$(handle(interface))" : ""
        interface_field = "$(typeof(interface.node)) $(interface.node.id) [$(findfirst(interface.node.interfaces, interface)):$(interface_handle)]"
        println(io, "$(string(entry_counter)): $(interface_field), $(string(msg_calc_func)) $(string(postproc))")
        entry_counter += 1
    end
end

function show(io::IO, nodes::Array{Node, 1})
     # Show node array (possibly an external schedule)
    println(io, "Nodes:")
    for entry in nodes
        println(io, "Node $(entry.id) of type $(typeof(entry))")
    end
end
