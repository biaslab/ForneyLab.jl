export GraphDump, pushNode!, pushEdge!, pushScore!, pushMessage!, pushMarginal!

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

function pushScore!(dump::GraphDump, id::Union{String, Symbol}, value::Float64, type::String)
    pushScore!(dump.steps[end], id, value, type)
    return dump
end

function pushMessage!(dump::GraphDump, edgeID::Union{String, Symbol}, type::String, message)
    pushMessage!(dump.steps[end], edgeID, type, message)
    return step
end

function pushMarginal!(dump::GraphDump, id::Union{String, Symbol}, edgeIDs::Vector{String}, marginal)
    pushMarginal!(dump.steps[end], id, edgeIDs, marginal)
    return step
end
