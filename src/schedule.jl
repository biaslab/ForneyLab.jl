export ScheduleEntry, Schedule, setPostProcessing!

type ScheduleEntry
    node::Node
    outbound_interface_id::Int64
    rule::Function  # Refers to the general message calculation rule; for example sumProduct! or vmp!.
    inbound_types::Vector{DataType}
    intermediate_outbound_type::DataType # Outbound type after update rule
    outbound_type::DataType # Outbound type after (optional) post-processing
    execute::Function # Compiled rule call: () -> rule(node, Val{outbound_interface_id}, rule_arguments...). Upon compilation execute() incorporates post-processing.
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

    isdefined(src, :inbound_types) && (duplicate.inbound_types = copy(src.inbound_types))
    isdefined(src, :intermediate_outbound_type) && (duplicate.outbound_type = src.intermediate_outbound_type)
    isdefined(src, :outbound_type) && (duplicate.outbound_type = src.outbound_type)
    isdefined(src, :post_processing) && (duplicate.post_processing = src.post_processing)

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

function generateScheduleByDFS!(outbound_interface::Interface,
                                backtrace::Vector{Interface} = Interface[],
                                call_list::Vector{Interface} = Interface[];
                                allowed_edges = false)
    # Private function to generate a sum product schedule by doing a DFS through the graph.
    # The graph is passed implicitly through the outbound_interface.
    #
    # IMPORTANT: the resulting schedule depends on the current messages stored in the factor graph.
    #
    # outbound_interface: find a schedule to calculate the outbound message on this interface
    # backtrace: backtrace for recursive implementation of DFS
    # call_list: holds the recursive calls
    # allowed_edges: either false or Set{Edge}. If a set is passed, the search will be restricted to edges in this set.
    #
    # Returns: Vector{Interface} (not an actual Schedule yet)

    node = outbound_interface.node

    # Apply stopping condition for recursion. When the same interface is called twice, this is indicative of an unbroken loop.
    if outbound_interface in call_list
        # Notify the user to break the loop with an initial message
        error("Loop detected around $(outbound_interface) Consider setting an initial message somewhere in this loop.")
    elseif outbound_interface in backtrace
        # This outbound_interface is already in the schedule
        return backtrace
    else # Stopping condition not satisfied
        push!(call_list, outbound_interface)
    end

    # Check all inbound messages on the other interfaces of the node
    for (id, interface) in enumerate(node.interfaces)
        is(interface, outbound_interface) && continue # Skip the outbound interface
        ( (typeof(allowed_edges)==Set{Edge}) && !(interface.edge in allowed_edges) ) && continue # Skip if the interface's edge is outside the allowed search set

        (interface.partner != nothing) || error("Disconnected interface should be connected: interface #$(interface_index) of $(typeof(node)) $(node.id)")

        if interface.partner.message == nothing # Required message missing, no breaker message is set
            if !(interface.partner in backtrace) # Don't recalculate stuff that's already in the schedule.
                generateScheduleByDFS!(interface.partner, backtrace, call_list, allowed_edges=allowed_edges) # Recursive call
            end
        end
    end

    # Update call_list and backtrace
    pop!(call_list)

    return push!(backtrace, outbound_interface)
end