"""
A `Cluster` specifies a collection of `edges` adjacent to `node` that belong to the same
`PosteriorFactor`. A joint marginal can be computed over a cluster.
"""
mutable struct Cluster <: Region
    id::Symbol
    node::FactorNode
    edges::Vector{Edge}

    function Cluster(node::FactorNode, edges::Vector{Edge})
        id = Symbol(join([edge.variable.id for edge in edges], "_"))
        self = new(id, node, edges)
        return self
    end
end

Base.isless(c1::Cluster, c2::Cluster) = isless("$(c1.id)", "$(c2.id)")
Base.isless(v::Variable, c::Cluster) = isless("$(v.id)", "$(c.id)")
Base.isless(c::Cluster, v::Variable) = isless("$(c.id)", "$(v.id)")

"""
Return the region that the node-edge combination belongs to (if available)
"""
function region(node::FactorNode, edge::Edge)
    dict = current_posterior_factorization.node_edge_to_cluster
    if haskey(dict, (node, edge))
        cl = dict[(node, edge)]
    else # No cluster is registered, return the edge variable
        cl = edge.variable
    end

    return cl::Region
end

"""
Return the local stochastic regions around `node`
"""
function localStochasticRegions(node::FactorNode)
    regions = Region[]
    for interface in node.interfaces
        partner = ultimatePartner(interface)
        if !isClamped(partner) # Exclude clamped
            push!(regions, region(node, interface.edge))
        end
    end

    return regions
end