export GraphData

struct GraphData
    nodes::Vector{NodeData}
    edges::Vector{EdgeData}
end

GraphData() = GraphData(Vector{NodeData}(), Vector{EdgeData}())

function pushNode!(data::GraphData, id, label, type, class, value=nothing)
    push!(data.nodes, NodeData(id, label, type, class, value))
    return data
end

function pushEdge!(data::GraphData, id, label, a, b, source, target)
    push!(data.edges, EdgeData(id, label, a, b, source, target))
    return data
end
