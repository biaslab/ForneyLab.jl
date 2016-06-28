export ScheduleEntry, Schedule, calculationRule

type ScheduleEntry{rule<:MessageCalculationRule}
    node::Node
    outbound_interface_id::Int64
    inbound_types::Vector{DataType} # Types of inbound messages/distributions.
    outbound_type::DataType         # Type of outbound distribution (the outbound distribution will be wrapped in a Message).
    approximation::DataType         # If the rule is an approximate one, this field specifies the type of approximation. Should be <: ApproximationType.
    unrolling_factor::Int64         # The number of factors that are unrolled.
    execute::Function               # Compiled rule call: () -> rule(node, Val{outbound_interface_id}, rule_arguments...). Invoked by execute(::ScheduleEntry).

    function ScheduleEntry(node::Node, outbound_interface_id::Int64)
        self = new{rule}(node, outbound_interface_id)
        self.unrolling_factor = 0 # can't leave this field undefined since its not a reference
        return self
    end

    function ScheduleEntry(interface::Interface)
        self = new{rule}(interface.node, findfirst(interface.node.interfaces, interface))
        self.unrolling_factor = 0 # can't leave this field undefined since its not a reference
        return self
    end
end

function Base.show{rule}(io::IO, entry::ScheduleEntry{rule})
    node = entry.node
    interface = node.interfaces[entry.outbound_interface_id]
    interface_handle = (handle(interface)!="") ? "($(handle(interface)))" : ""
    approx = isdefined(entry, :approximation) ? "(Approx.: $(entry.approximation)) " : ""
    println(io, replace("$(approx)$(rule) on $(typeof(node)) $(interface.node.id) interface $(entry.outbound_interface_id) $(interface_handle)", "ForneyLab.", ""))
    if isdefined(entry, :inbound_types) && isdefined(entry, :outbound_type)
        println(io, replace("$(entry.inbound_types) -> Message{$(entry.outbound_type)}", "ForneyLab.", ""))
    end
    if entry.unrolling_factor > 0
        println(io, "The partitioned distributions are unrolled and the rule is applied to each of the $(entry.unrolling_factor) factors.")
    end
end

calculationRule{rule}(::ScheduleEntry{rule}) = rule

Base.deepcopy(::ScheduleEntry) = error("deepcopy(::ScheduleEntry) is not possible. You should construct a new ScheduleEntry or use copy(::ScheduleEntry).")

function Base.copy(src::ScheduleEntry)
    duplicate = typeof(src)(src.node, src.outbound_interface_id)

    duplicate.unrolling_factor = src.unrolling_factor
    isdefined(src, :inbound_types) && (duplicate.inbound_types = copy(src.inbound_types))
    isdefined(src, :outbound_type) && (duplicate.outbound_type = src.outbound_type)
    isdefined(src, :approximation) && (duplicate.approximation = src.approximation)

    return duplicate
end

function setOutboundType!(entry::ScheduleEntry, outbound_type::DataType)
    if outbound_type <: ProbabilityDistribution
        entry.outbound_type = outbound_type
    elseif outbound_type <: Approximation
        entry.outbound_type = outbound_type.parameters[1]
        entry.approximation = outbound_type.parameters[2]
    else
        error("Invalid message type specification: $(outbound_type). Should be subtype of ProbabilityDistribution or Approximation.")
    end

    return entry
end

function convert{rule<:MessageCalculationRule}(::Type{ScheduleEntry{rule}}, interface::Interface)
    node = interface.node
    interface_id = findfirst(node.interfaces, interface)
    return ScheduleEntry{rule}(node, interface_id)
end

"""
Defines a sequence of Message calculations
"""
typealias Schedule Vector{ScheduleEntry}

Base.deepcopy{T<:ScheduleEntry}(src::Vector{T}) = ScheduleEntry[copy(entry) for entry in src]

function convert{rule<:MessageCalculationRule}(::Type{Schedule}, interfaces::Vector{Interface}, ::Type{rule})
    return ScheduleEntry[convert(ScheduleEntry{rule}, interface) for interface in interfaces]
end

function Base.show(io::IO, schedule::Vector{ScheduleEntry})
    println(io, "Message passing schedule")
    println(io, "------------------------\n")
    for i=1:length(schedule)
        println("$(i).")
        show(schedule[i])
        println("")
    end
end

# TODO: remove
function generateScheduleByDFS!(outbound_interface::Interface,
                                backtrace::Vector{Interface} = Interface[],
                                call_list::Vector{Interface} = Interface[];
                                allowed_edges::Set{Edge}      = Set{Edge}(),
                                breaker_sites::Set{Interface} = Set{Interface}())

    # Internal function to generate a sum product schedule by doing a DFS through the graph.
    # The graph is passed implicitly through the outbound_interface.
    #
    # outbound_interface:   find a schedule to calculate the outbound message on this interface
    # backtrace:            backtrace for recursive implementation of DFS
    # call_list:            holds the recursive calls
    # allowed_edges:        a non-empty set indicates that the search will be restricted to the edges from this set
    # breaker_sites:        terminate the DFS at these sites
    #
    # Returns: Vector{Interface}

    node = outbound_interface.node

    # Apply stopping condition for recursion. When the same interface is called twice, this is indicative of an unbroken loop.
    if outbound_interface in call_list
        # Notify the user to break the loop with a breaker message
        error("Loop detected around $(outbound_interface). Consider using a loopy algorithm and setting a breaker message somewhere in this loop.")
    elseif outbound_interface in backtrace
        # This outbound_interface is already in the schedule
        return backtrace
    else # Stopping condition not satisfied
        push!(call_list, outbound_interface)
    end

    # Check all inbound messages on the other interfaces of the node
    for (id, interface) in enumerate(node.interfaces)
        is(interface, outbound_interface) && continue # Skip the outbound interface
        ( (length(allowed_edges) > 0) && !(interface.edge in allowed_edges) ) && continue # Skip if the interface's edge is outside the allowed search set

        (interface.partner != nothing) || error("Disconnected interface should be connected: $(interface)")

        if !(interface.partner in breaker_sites) # No termination because of breaker message
            if !(interface.partner in backtrace) # Don't recalculate stuff that's already in the schedule.
                generateScheduleByDFS!(interface.partner, backtrace, call_list, allowed_edges=allowed_edges, breaker_sites=breaker_sites) # Recursive call
            end
        end
    end

    # Update call_list and backtrace
    pop!(call_list)

    return push!(backtrace, outbound_interface)
end
