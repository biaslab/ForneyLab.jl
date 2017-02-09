export FactorNode

abstract FactorNode
abstract DeltaFactor <: FactorNode
abstract SoftFactor <: FactorNode

Base.isless(n1::FactorNode, n2::FactorNode) = isless("$(n1.id)", "$(n2.id)")

show(io::IO, node::FactorNode) = println(io, "$(typeof(node)) with id $(node.id)")

function show(io::IO, nodes::Union{Vector{FactorNode},Set{FactorNode}})
    println(io, "FactorNodes:")
    for node in nodes
        show(io, node)
    end
end