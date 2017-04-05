export
MessageUpdateRule,
ScheduleEntry,
Schedule

"""
A `ScheduleEntry` defines a message computation.
The `msg_update_rule <: MessageUpdateRule` defines the rule that is used
to calculate the message coming out of `interface`.
"""
type ScheduleEntry
    interface::Interface
    msg_update_rule::DataType
end


typealias Schedule Vector{ScheduleEntry}


"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageUpdateRule


"""
`summaryPropagationSchedule(variables)` builds a generic summary propagation
`Schedule` for calculating the marginal distributions of every variable in
`variables`. The message update rule in the schedule entries is set to `Void`.
"""
function summaryPropagationSchedule(variables::Vector{Variable})
    # We require the marginal distribution of every variable in variables.
    # If a variable relates to multiple edges, this indicates an equality constraint.
    # Therefore, we only need to consider one arbitrary edge to calculate the marginal.
    # For the sake of consistency, we always take the first edge.
    seed_interfaces = Interface[]
    for variable in variables
        edge = first(variable.edges)
        (edge.a != nothing) && push!(seed_interfaces, edge.a)
        (edge.b != nothing) && push!(seed_interfaces, edge.b)
    end

    # Determine a feasible ordering of message updates
    dg = summaryDependencyGraph(current_graph)
    iface_list = children(unique(seed_interfaces), dg)
    # Build a schedule; Void indicates an unspecified message update rule
    schedule = [ScheduleEntry(iface, Void) for iface in iface_list]

    return schedule
end
