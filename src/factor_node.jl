export FactorNode

abstract type FactorNode end
abstract type DeltaFactor <: FactorNode end
abstract type SoftFactor <: FactorNode end

Base.isless(n1::FactorNode, n2::FactorNode) = isless("$(n1.id)", "$(n2.id)")

show(io::IO, node::FactorNode) = println(io, "$(typeof(node)) with id $(node.id)")

function show(io::IO, nodes::Union{Vector{FactorNode},Set{FactorNode}})
    println(io, "FactorNodes:")
    for node in nodes
        show(io, node)
    end
end

slug(node_type::Type{T}) where T<:FactorNode = replace(string(T), "ForneyLab." => "")