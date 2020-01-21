export
MarginalUpdateRule,
MarginalEntry,
MarginalTable,
marginalTable

"""
A `MarginalUpdateRule` specifies how a (joint) marginal is calculated from
incoming messages (and a node function).
"""
abstract type MarginalUpdateRule end
abstract type Product <: MarginalUpdateRule end

"""
A `MarginalEntry` defines a marginal computation.
The `marginal_update_rule <: MarginalUpdateRule` defines the rule that is used
to calculate the (joint) marginal over `target`.
"""
mutable struct MarginalEntry
    target::Union{Variable, Cluster}
    interfaces::Vector{Interface}
    marginal_update_rule::DataType

    # Fields for algorithm assembly
    marginal_id::Symbol # Specify the marginal identifier
    inbounds::Vector{Any} # Specify the inbounds required for the marginal update

    MarginalEntry() = new()
    MarginalEntry(target::Union{Variable, Cluster}, interfaces::Vector{Interface}, marginal_update_rule::DataType) = new(target, interfaces, marginal_update_rule)
end

"""
A `MarginalTable` defines the update order for marginal computations.
"""
const MarginalTable = Vector{MarginalEntry}

"""
Generate a `MarginalTable` that computes the marginals for target variables
"""
function marginalTable(targets::Vector{Variable})
    marginal_table = MarginalEntry[]
    for target in targets
        edge = first(target.edges) # For the sake of consistency, we always take the first edge as reference point for marginal computations.
        if edge.a == nothing # First handle cases where there is a `dangling` edge
            push!(marginal_table, MarginalEntry(target, [edge.b], Nothing))
        elseif edge.b == nothing
            push!(marginal_table, MarginalEntry(target, [edge.a], Nothing))
        elseif isa(edge.a.node, Clamp) || isa(edge.b.node, Clamp)
            continue # Do not compute marginals for clamped edges
        else
            push!(marginal_table, MarginalEntry(target, [edge.a, edge.b], Product))
        end
    end

    return marginal_table
end

marginalTable(target::Variable) = marginalTable([target])

"""
Generate a mapping from target to marginal entry.
"""
function targetToMarginalEntry(table::MarginalTable)
    mapping = Dict{Union{Cluster, Variable}, MarginalEntry}()
    for entry in table
        target = entry.target
        mapping[target] = entry
    end

    return mapping
end