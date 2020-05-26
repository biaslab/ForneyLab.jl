export
Message,
MessageUpdateRule,
ScheduleEntry,
Schedule

"""Encodes a message, which is a probability distribution with a scaling factor"""
struct Message{family<:FactorFunction, var_type<:VariateType} # Note that parameter order is switched w.r.t. ProbabilityDistribution, for ease of overloading
    dist::ProbabilityDistribution{var_type, family}
    scaling_factor::Any

    Message{F, V}(dist::ProbabilityDistribution{V, F}) where {F, V}= new(dist) # Constructor for unspecified scaling factor
end

Message(dist::ProbabilityDistribution{V, F}) where {F<:FactorFunction, V<:VariateType} = Message{F, V}(dist)

Message(var_type::Type{V}, family::Type{F}; kwargs...) where {F<:FactorFunction, V<:VariateType} = Message{family, var_type}(ProbabilityDistribution(var_type, family; kwargs...))

function Message(family::Type{F}; kwargs...) where F
    dist = ProbabilityDistribution(family; kwargs...)
    var_type = variateType(dist)

    return Message{family, var_type}(dist)
end

family(msg_type::Type{Message{F}}) where F<:FactorFunction = F

function show(io::IO, msg::Message)
    if isdefined(msg, :scaling_factor)
        println(io, "Message: 1/$(format(msg.scaling_factor)) * $(format(msg.dist))")
    else
        println(io, "Message: $(format(msg.dist))")
    end
end

"""Special inheritance rules for parametric Message types"""
matches(::Type{T}, ::Type{T}) where T<:Message = true
matches(::Type{Message{F, V}}, ::Type{Message}) where {F<:FactorFunction, V<:VariateType} = true
matches(::Type{Message{F}}, ::Type{Message}) where {F<:FactorFunction} = true
matches(Ta::Type{Message{Fa, Va}}, Tb::Type{Message{Fb, Vb}}) where {Fa<:FactorFunction, Fb<:FactorFunction, Va<:VariateType, Vb<:VariateType} = (Va==Vb) && (Fa<:Fb)
matches(Ta::Type{Message{Fa, Va}}, Tb::Type{Message{Fb}}) where {Fa<:FactorFunction, Fb<:FactorFunction, Va<:VariateType} = (Fa<:Fb)
matches(Ta::Type{Message{Fa}}, Tb::Type{Message{Fb}}) where {Fa<:FactorFunction, Fb<:FactorFunction} = (Fa<:Fb)
matches(::Type{Nothing}, ::Type{T}) where T<:Message = false
matches(::Type{P}, ::Type{M}) where {P<:ProbabilityDistribution, M<:Message} = false
matches(::Type{M}, ::Type{P}) where {P<:ProbabilityDistribution, M<:Message} = false

function ==(t::Message{fam_t, var_t}, u::Message{fam_u, var_u}) where {fam_t<:FactorFunction, var_t<:VariateType, fam_u<:FactorFunction, var_u<:VariateType}
    (fam_t == fam_u) || return false
    (var_t == var_u) || return false
    (t.dist == u.dist) || return false
    if isdefined(t, :scaling_factor) && isdefined(u, :scaling_factor)
        (t.scaling_factor == u.scaling_factor) || return false
    end
    isdefined(t, :scaling_factor) && !isdefined(u, :scaling_factor) && return false
    !isdefined(t, :scaling_factor) && isdefined(u, :scaling_factor) && return false

    return true
end

"""
A MessageUpdateRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageUpdateRule)` to list the available rules.
"""
abstract type MessageUpdateRule end

"""
A `ScheduleEntry` defines a message computation.
The `message_update_rule <: MessageUpdateRule` defines the rule that is used
to calculate the message coming out of `interface`.
"""
mutable struct ScheduleEntry
    interface::Interface
    message_update_rule::Type
    internal_schedule::Vector{ScheduleEntry}

    # Fields for algorithm assembly
    schedule_index::Int # Specify the position in the schedule
    initialize::Bool # Indicate whether this mesage is initialized before schedule execution
    family::Type # Specify the outbound family for this entry
    dimensionality::Tuple # Specify the outbound dimensionality
    inbounds::Vector{Any} # Specify the inbounds required for computing the message update

    ScheduleEntry() = new()
    ScheduleEntry(interface::Interface, message_update_rule::Type) = new(interface, message_update_rule)
end

function show(io::IO, entry::ScheduleEntry)
    rule_str = replace(string(entry.message_update_rule), "ForneyLab." => "") # Remove "Forneylab."
    internal_schedule = isdefined(entry, :internal_schedule) ? "(INTERNAL SCHEDULE) " : ""
    print(io, "$(internal_schedule)$(rule_str) on $(entry.interface)")
end

function ==(a::ScheduleEntry, b::ScheduleEntry)
    return (a.interface == b.interface) && (a.message_update_rule == b.message_update_rule)
end

"""
A `Schedule` defines the update order for message computations.
"""
const Schedule = Vector{ScheduleEntry}

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
                (node_interface === interface.partner) && continue
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
`summaryPropagationSchedule(variables, clusters)` builds a generic
`Schedule` for calculating the marginal distributions of every variable in
`variables` and every cluster in `clusters`. The message update rule in each 
schedule entry is set to `Nothing`.
"""
function summaryPropagationSchedule(variables::Vector{Variable}, clusters::Vector{Cluster}=Cluster[]; limit_set=edges(current_graph), target_sites=Interface[], breaker_sites=Interface[])
    # We require the marginal distribution of every variable in variables.
    # If a variable relates to multiple edges, this indicates an equality constraint.
    # Therefore, we only need to consider one arbitrary edge to calculate the marginal.
    for variable in variables
        edge = first(variable.edges) # For the sake of consistency, we always take the first edge.
        (edge.a != nothing) && !isa(edge.a.node, Terminal) && !isa(edge.a.node, Clamp) && push!(target_sites, edge.a)
        (edge.b != nothing) && !isa(edge.b.node, Terminal) && !isa(edge.b.node, Clamp) && push!(target_sites, edge.b)
    end

    # A cluster indicates the computation of a joint marginal belief.
    # In this case, all incoming messages on internal edges are required.
    for cluster in clusters
        for iface in cluster.node.interfaces
            partner = ultimatePartner(iface)
            (iface.edge in limit_set) && (partner != nothing) && !isa(partner.node, Terminal) && !isa(partner.node, Clamp) && push!(target_sites, partner)
        end
    end

    # Determine a feasible ordering of message updates
    dg = summaryDependencyGraph(limit_set)
    iface_list = children(unique(target_sites), dg, breaker_sites=Set(breaker_sites))
    # Build a schedule; Nothing indicates an unspecified message update rule
    schedule = [ScheduleEntry(iface, Nothing) for iface in iface_list]

    return schedule
end

"""
inferUpdateRules!(schedule) infers specific message update rules for all schedule entries.
"""
function inferUpdateRules!(schedule::Schedule; inferred_outbound_types=Dict{Interface, Type}())
    for entry in schedule
        (entry.message_update_rule == Nothing) && error("No msg update rule type specified for $(entry)")
        if !isconcretetype(entry.message_update_rule)
            # In this case entry.message_update_rule is a update rule type, but not a specific rule.
            # Here we infer the specific rule that is applicable, which should be a subtype of entry.message_update_rule.
            inferUpdateRule!(entry, entry.message_update_rule, inferred_outbound_types)
        end
        # Store the rule's outbound type
        inferred_outbound_types[entry.interface] = outboundType(entry.message_update_rule)
    end

    return schedule
end


"""
Flatten a schedule by inlining all internal message passing schedules.
This yields a simple, linear schedule.
"""
function flatten(schedule::Schedule)
    # Check if there are any internal schedules to inline
    if !any([isdefined(entry, :internal_schedule) for entry in schedule])
        return schedule
    end

    flat_schedule = ScheduleEntry[]
    passed_interfaces = Dict{Interface,Bool}()
    for entry in schedule
        if isdefined(entry, :internal_schedule)
            flat_subschedule = flatten(entry.internal_schedule)
            for internal_entry in flat_subschedule
                if !haskey(passed_interfaces, internal_entry.interface)
                    push!(flat_schedule, internal_entry)
                    passed_interfaces[internal_entry.interface] = true
                end
            end
        else
            if !haskey(passed_interfaces, entry.interface)
                push!(flat_schedule, entry)
                passed_interfaces[entry.interface] = true
            end
        end
    end

    return flat_schedule
end


"""
Contruct a condensed schedule.
"""
function condense(schedule::Schedule)
    condensed_schedule = ScheduleEntry[]
    for entry in schedule
        if isdefined(entry, :internal_schedule)
            entry.internal_schedule = condense(entry.internal_schedule)
        elseif !isa(entry.interface.node, Clamp)
            push!(condensed_schedule, entry)
        end
    end

    return condensed_schedule
end


"""
Generate a mapping from interface to schedule entry.
Multiple interfaces can map to the same schedule entry if the
graph contains composite nodes.
"""
function interfaceToScheduleEntry(schedule::Schedule)
    mapping = Dict{Interface, ScheduleEntry}()
    for entry in schedule
        interface = entry.interface
        mapping[interface] = entry
        while (interface.partner != nothing) && isa(interface.partner.node, Terminal)
            interface = interface.partner.node.outer_interface
            mapping[interface] = entry
        end
    end

    return mapping
end
