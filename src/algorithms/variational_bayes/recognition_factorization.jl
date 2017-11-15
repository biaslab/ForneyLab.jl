import Base.factor
export factor, RecognitionFactor, RecognitionFactorization, currentRecognitionFactorization

"""
A RecognitionFactor specifies the subset of variables that comprise
a joint factor in the recognition factorization.
"""
type RecognitionFactor
    id::Symbol
    variables::Set{Variable}
    internal_edges::Set{Edge}

    function RecognitionFactor(variables::Set{Variable}; rfz=currentRecognitionFactorization(), id=generateId(RecognitionFactor))
        # Collect variables on internal edges connected to external nodes
        internal_edges = extend(edges(variables))
        internal_edges_connected_to_external_nodes = intersect(edges(nodes(internal_edges)), internal_edges)
        recognition_variables = Set{Variable}([edge.variable for edge in internal_edges_connected_to_external_nodes])

        self = new(id, union(variables, recognition_variables), internal_edges)
        rfz.recognition_factors[id] = self # Register new factor with recognition factorization

        return self 
    end
end

RecognitionFactor(variable::Variable; id=generateId(RecognitionFactor)) = RecognitionFactor(Set([variable]), id=id)
RecognitionFactor(variables::Vector{Variable}; id=generateId(RecognitionFactor)) = RecognitionFactor(Set(variables), id=id)

function draw(rf::RecognitionFactor; schedule=ScheduleEntry[], args...)
    subgraph_nodes = nodes(rf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), rf.internal_edges)
    ForneyLab.graphviz(ForneyLab.genDot(subgraph_nodes, rf.internal_edges, schedule=schedule, external_edges=external_edges); args...)
end

"""
Find the smallest legal subgraph (connected through deterministic nodes) that includes the argument edges
"""
function extend(edge_set::Set{Edge})
    cluster = Set{Edge}() # Set to fill with edges in cluster
    edges = copy(edge_set)
    while length(edges) > 0 # As long as there are unchecked edges connected through deterministic nodes
        current_edge = pop!(edges) # Pick one
        push!(cluster, current_edge) # Add to edge cluster
        for node in [current_edge.a.node, current_edge.b.node] # Check both head and tail node for deterministic type
            if isa(node, DeltaFactor)
                for interface in node.interfaces
                    if !is(interface.edge, current_edge) && !(interface.edge in cluster) # Is next level edge not seen yet?
                        push!(edges, interface.edge) # Add to buffer to visit sometime in the future
                    end
                end
            end
        end
    end

    return cluster
end

type RecognitionFactorization
    graph::FactorGraph
    recognition_factors::Dict{Symbol, RecognitionFactor}
end

"""
Return currently active RecognitionFactorization.
Create one if there is none.
"""
function currentRecognitionFactorization()
    try
        return current_recognition_factorization
    catch
        return RecognitionFactorization()
    end
end

setCurrentRecognitionFactorization(rf::RecognitionFactorization) = global current_recognition_factorization = rf

RecognitionFactorization() = setCurrentRecognitionFactorization(RecognitionFactorization(currentGraph(), Dict{Symbol, RecognitionFactor}()))

function nodesConnectedToExternalEdges(recognition_factor::RecognitionFactor)
    internal_edges = recognition_factor.internal_edges
    subgraph_nodes = nodes(internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)

    return nodes_connected_to_external_edges
end