export GraphDump

mutable struct GraphDump
    graph::GraphData
    marginals::Vector{AlgorithmStep}
end

GraphDump() = GraphDump(GraphData(Vector{NodeData}(),
                                  Vector{EdgeData}()),
                        Vector{AlgorithmStep}())
