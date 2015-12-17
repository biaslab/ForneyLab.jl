export ScheduleEntry, Schedule, setPostProcessing!

type ScheduleEntry
    node::Node
    outbound_interface_id::Int64
    rule::Function  # Refers to the general message calculation rule; for example sumProduct! or vmp!.
    inbound_types::Vector{DataType}
    outbound_type::Type{ProbabilityDistribution}
    execute::Function # Compiled rule call: () -> rule(node, Val{outbound_interface_id}, rule_arguments...)
    post_processing::Function

    function ScheduleEntry(node::Node, outbound_interface_id::Int64, rule::Function)
        return self = new(node, outbound_interface_id, rule)
    end
end

function ScheduleEntry(node::Node, outbound_interface_id::Int64, rule::Function, post_processing::Function)
    schedule_entry = ScheduleEntry(node, outbound_interface_id, rule)
    schedule_entry.post_processing = post_processing

    return schedule_entry
end

Base.deepcopy(::ScheduleEntry) = error("deepcopy(::ScheduleEntry) is not possible. You should construct a new ScheduleEntry or use copy(::ScheduleEntry).")

function Base.copy(src::ScheduleEntry)
    duplicate = ScheduleEntry(src.node, src.outbound_interface_id, src.rule)
    if isdefined(src, :inbound_types)
        duplicate.inbound_types = copy(src.inbound_types)
    end
    if isdefined(src, :outbound_type)
        duplicate.outbound_type = src.outbound_type
    end
    if isdefined(src, :post_processing)
        duplicate.post_processing = src.post_processing
    end

    return duplicate
end

function setPostProcessing!(schedule_entry::ScheduleEntry, post_processing::Function)
    schedule_entry.post_processing = post_processing
    return schedule_entry
end

typealias Schedule Array{ScheduleEntry, 1}

Base.deepcopy(src::Schedule) = ScheduleEntry[copy(entry) for entry in src]

function setPostProcessing!(schedule::Schedule, interface::Interface, post_processing::Function)
    for entry in schedule
        if is(entry.node.interfaces[entry.outbound_interface_id], interface)
            entry.post_processing = post_processing # Edit in place
            return schedule
        end
    end
end

# Convert interfaces to schedule
function convert(::Type{ScheduleEntry}, interface::Interface, rule::Function = sumProduct!)
    node = interface.node
    interface_id = findfirst(node.interfaces, interface)
    return ScheduleEntry(node, interface_id, rule)
end

function convert(::Type{Schedule}, interfaces::Vector{Interface}, rule::Function)
    return ScheduleEntry[convert(ScheduleEntry, interface, rule) for interface in interfaces]
end

function show(io::IO, schedule_entry::ScheduleEntry)
    node = schedule_entry.node
    interface = node.interfaces[schedule_entry.outbound_interface_id]
    interface_handle = (handle(interface)!="") ? "($(handle(interface)))" : ""
    println(io, replace("$(schedule_entry.rule) on $(typeof(node)) $(interface.node.id) interface $(schedule_entry.outbound_interface_id) $(interface_handle)", "ForneyLab.", ""))
    if isdefined(schedule_entry, :inbound_types) && isdefined(schedule_entry, :outbound_type)
        println(io, replace("$(schedule_entry.inbound_types) -> Message{$(schedule_entry.outbound_type)}", "ForneyLab.", ""))
    end
    if isdefined(schedule_entry, :post_processing)
        println(io, replace("Post processing: $(schedule_entry.post_processing)", "ForneyLab.", ""))
    end
end

function show(io::IO, schedule::Schedule)
    println(io, "Message passing schedule")
    println(io, "-----------------------------------------------")
    for i=1:length(schedule)
        println("$(i).")
        show(schedule[i])
        println("")
    end
end