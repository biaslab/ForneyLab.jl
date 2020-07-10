export GraphDump

struct GraphDump
    graph::GraphData
    steps::Vector{AlgorithmStep}
end

function GraphDump()
    d = GraphDump(GraphData(), Vector{AlgorithmStep}())
    push!(d.steps, AlgorithmStep())
    return d
end

function pushNode!(dump::GraphDump, id, label, type, class, value=nothing)
    pushNode!(dump.graph, id, label, type, class, value)
    return dump
end

function pushEdge!(dump::GraphDump, id, label, a, b, source, target)
    pushEdge!(dump.graph, id, label, a, b, source, target)
    return dump
end
