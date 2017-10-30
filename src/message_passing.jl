export
Message,
MessageUpdateRule,
ScheduleEntry,
Schedule

import Base: ==

"""Encodes a message, which is a probability distribution with a scaling factor"""
immutable Message{T<:ProbabilityDistribution}
    dist::T
    scaling_factor::Any

    Message{U<:ProbabilityDistribution}(dist::U) = new{U}(dist) # Constructor for unspecified scaling factor
end

Message{T}(dist::T) = Message{T}(dist)

function =={T<:ProbabilityDistribution, U<:ProbabilityDistribution}(t::Message{T}, u::Message{U})
    (T == U) || return false
    (t.dist == u.dist) || return false
    if isdefined(t, :scaling_factor) && isdefined(u, :scaling_factor)
        (t.scaling_factor == u.scaling_factor) || return false
    end
    isdefined(t, :scaling_factor) && !isdefined(u, :scaling_factor) && return false
    !isdefined(t, :scaling_factor) && isdefined(u, :scaling_factor) && return false
    
    return true
end

"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageUpdateRule

"""
A `ScheduleEntry` defines a message computation.
The `msg_update_rule <: MessageUpdateRule` defines the rule that is used
to calculate the message coming out of `interface`.
"""
type ScheduleEntry
    interface::Interface
    msg_update_rule::DataType
end

function show(io::IO, entry::ScheduleEntry)
    rule_str = replace(string(entry.msg_update_rule), "ForneyLab.", "") # Remove "Forneylab."
    print(io, "$(rule_str) on $(entry.interface)")
end

function ==(a::ScheduleEntry, b::ScheduleEntry)
    return (a.interface == b.interface) && (a.msg_update_rule == b.msg_update_rule)
end

typealias Schedule Vector{ScheduleEntry}

function show(io::IO, schedule::Schedule)
    condensed_schedule = condense(schedule)
    idx = 1
    for entry in schedule
        if entry in condensed_schedule
            print(io, "$idx.\t$(entry)")
            idx += 1
        else
            print(io, "\t$(entry)")
        end
    end
end

"""
summaryDependencyGraph(edgeset)

Returns a DependencyGraph (directed graph) that encodes the dependencies among
summary messages (such as sum-product messages) in `edgeset`.
All Interfaces in `edgeset` are vertices in the dependency graph.
The dependency graph can be used for loop detection, scheduling, etc.
"""
function summaryDependencyGraph(edgeset::Set{Edge}; reverse_edges=false)
    # Create dependency graph object
    dg = DependencyGraph{Interface}()

    # Add all Interfaces in edgeset as vertices in dg
    for edge in sort(collect(edgeset))
        isa(edge.a, Interface) && addVertex!(dg, edge.a)
        isa(edge.b, Interface) && addVertex!(dg, edge.b)
    end

    # Add all summary dependencies
    for interface in dg.vertices
        if isa(interface.partner, Interface) # interface is connected to an Edge
            for node_interface in interface.partner.node.interfaces
                is(node_interface, interface.partner) && continue
                (node_interface.edge in edgeset) || continue
                if reverse_edges
                    addEdge!(dg, interface, node_interface)
                else
                    addEdge!(dg, node_interface, interface)
                end
            end
        end
    end

    return dg
end

"""
`summaryPropagationSchedule(variables)` builds a generic summary propagation
`Schedule` for calculating the marginal distributions of every variable in
`variables`. The message update rule in each schedule entry is set to `Void`.
"""
function summaryPropagationSchedule(variables::Vector{Variable}; limit_set=edges(current_graph), target_sites=Interface[], breaker_sites=Interface[])
    # We require the marginal distribution of every variable in variables.
    # If a variable relates to multiple edges, this indicates an equality constraint.
    # Therefore, we only need to consider one arbitrary edge to calculate the marginal.
    for variable in variables
        edge = first(variable.edges) # For the sake of consistency, we always take the first edge.
        (edge.a != nothing) && push!(target_sites, edge.a)
        (edge.b != nothing) && push!(target_sites, edge.b)
    end

    # Determine a feasible ordering of message updates
    dg = summaryDependencyGraph(limit_set)
    iface_list = children(unique(target_sites), dg, breaker_sites=Set(breaker_sites))
    # Build a schedule; Void indicates an unspecified message update rule
    schedule = [ScheduleEntry(iface, Void) for iface in iface_list]
    
    return schedule
end

summaryPropagationSchedule(variable::Variable; limit_set=edges(current_graph), target_sites=Interface[], breaker_sites=Interface[]) = summaryPropagationSchedule([variable], limit_set=limit_set, target_sites=target_sites, breaker_sites=breaker_sites)

"""
inferUpdateRules!(schedule) infers specific message update rules for all schedule entries.
"""
function inferUpdateRules!(schedule::Schedule; inferred_outbound_types=Dict{Interface, DataType}())
    # Dict to hold all inferred message types
    for entry in schedule
        (entry.msg_update_rule == Void) && error("No msg update rule type specified for $(entry)")
        if !isleaftype(entry.msg_update_rule)
            # In this case entry.msg_update_rule is a update rule type, but not a specific rule.
            # Here we infer the specific rule that is applicable, which should be a subtype of entry.msg_update_rule.
            inferUpdateRule!(entry, entry.msg_update_rule, inferred_outbound_types)
        end
        # Store the rule's outbound type
        inferred_outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    return schedule
end

"""
Contruct a condensed schedule.
"""
function condense(schedule::Schedule)
    condensed_schedule = ScheduleEntry[]
    for entry in schedule
        if !isa(entry.interface.node, Clamp)
            push!(condensed_schedule, entry)
        end
    end

    return condensed_schedule
end
