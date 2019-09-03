export
MarginalUpdateRule,
MarginalScheduleEntry,
MarginalSchedule,
marginalSchedule

abstract type AbstractCluster end

"""
A `MarginalUpdateRule` specifies how a (joint) marginal is calculated from
incoming messages (and a node function).
"""
abstract type MarginalUpdateRule end
abstract type Product <: MarginalUpdateRule end

"""
A `MarginalScheduleEntry` defines a marginal computation.
The `marginal_update_rule <: MarginalUpdateRule` defines the rule that is used
to calculate the (joint) marginal over `target`.
"""
mutable struct MarginalScheduleEntry
    target::Union{Variable, AbstractCluster}
    interfaces::Vector{Interface}
    marginal_update_rule::DataType
end

"""
A `MarginalSchedule` defines the update order for marginal computations.
"""
const MarginalSchedule = Vector{MarginalScheduleEntry}

"""
Generate a `MarginalSchedule` that computes the marginals for target variables
"""
function marginalSchedule(targets::Vector{Variable})
    marginal_schedule = MarginalScheduleEntry[]
    for target in targets
        edge = first(target.edges) # For the sake of consistency, we always take the first edge as reference point for marginal computations.
        if edge.a == nothing # First handle cases where there is a `dangling` edge
            push!(marginal_schedule, MarginalScheduleEntry(target, [edge.b], Nothing))
        elseif edge.b == nothing
            push!(marginal_schedule, MarginalScheduleEntry(target, [edge.a], Nothing))
        elseif isa(edge.a.node, Clamp) || isa(edge.b.node, Clamp)
            continue # Do not compute marginals for clamped edges
        else
            push!(marginal_schedule, MarginalScheduleEntry(target, [edge.a, edge.b], Product))
        end
    end

    return marginal_schedule
end

marginalSchedule(target::Variable) = marginalSchedule([target])