export
MarginalUpdateRule,
MarginalScheduleEntry,
MarginalSchedule

abstract AbstractCluster

"""
A MarginalUpdateRule specifies how a (joint) marginal is calculated from
incoming messages (and a node function).
"""
abstract MarginalUpdateRule
abstract MarginalRule{factor_type} <: MarginalUpdateRule
abstract Product <: MarginalUpdateRule

"""
A `MarginalScheduleEntry` defines a marginal computation.
The `marginal_update_rule <: MarginalUpdateRule` defines the rule that is used
to calculate the (joint) marginal over `target`.
"""
type MarginalScheduleEntry
    target::Union{Variable, AbstractCluster}
    interfaces::Vector{Interface}
    marginal_update_rule::DataType
end

typealias MarginalSchedule Vector{MarginalScheduleEntry}

"""
Construct a MarginalScheduleEntry for computing the marginal over `variable`
through multiplication of colliding messages.
"""
function MarginalScheduleEntry(variable::Variable)
    edge = first(variable.edges) # For the sake of consistency, we always take the first edge as reference point for marginal computations.
    if edge.a == nothing # First handle cases where there is a `dangling` edge
        entry = MarginalScheduleEntry(variable, [edge.b], Void)
    elseif edge.b == nothing
        entry = MarginalScheduleEntry(variable, [edge.a], Void)
    else
        entry = MarginalScheduleEntry(variable, [edge.a, edge.b], Product)
    end

    return entry
end

"""
marginalSchedule() generates a marginal schedule that computes the marginals for each target entry
"""
marginalSchedule(targets::Vector{Variable}) = [MarginalScheduleEntry(target) for target in targets]