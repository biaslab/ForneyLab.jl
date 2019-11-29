"""
A `Cluster` specifies a collection of `edges` adjacent to `node` that belong to the same
`RecognitionFactor`. A joint marginal can be computed over a cluster.
"""
mutable struct Cluster <: AbstractCluster
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

"""
Return the cluster that the node-edge combination belongs to (if available)
"""
function cluster(node::FactorNode, edge::Edge)
    dict = current_algorithm.node_edge_to_cluster
    if haskey(dict, (node, edge))
        cl = dict[(node, edge)]
    else # No cluster is found, return the variable itself
        cl = edge.variable
    end

    return cl::Union{Cluster, Variable}
end

"""
Return the local clusters/variables around `node`
"""
localClusters(node::FactorNode) = [cluster(node, interface.edge) for interface in node.interfaces]