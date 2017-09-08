export RecognitionFactor

"""
A RecognitionFactor specifies the subset of variables that comprise
a joint factor in the recognition factorization.
"""
type RecognitionFactor
    id::Symbol
    internal_edges::Set{Edge}
    recognition_distributions::Dict{Variable, ProbabilityDistribution}

    function RecognitionFactor(variables::Set{Variable}, dist::ProbabilityDistribution; id=generateId(RecognitionFactor))
        internal_edges = extend(edges(variables))
        recognition_distributions = Dict{Variable, ProbabilityDistribution}()
        for variable in variables
            recognition_distributions[variable] = dist
        end

        return new(id, internal_edges, recognition_distributions)
    end
end

RecognitionFactor(variable, dist; id=generateId(RecognitionFactor)) = RecognitionFactor(Set([variable]), dist, id=id)

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

# type RecognitionFactorization
#     recognition_factors::Vector{RecognitionFactor}
# end

# """
# Return currently active RecognitionFactorization.
# Create one if there is none.
# """
# function currentRecognitionFactorization()
#     try
#         return current_recognition_factorization
#     catch
#         return RecognitionFactorization()
#     end
# end

# setCurrentRecognitionFactorization(rf::RecognitionFactorization) = global current_recognition_factorization = rf

# RecognitionFactorization() = setCurrentRecognitionFactorization(RecognitionFactorization(RecognitionFactor[]))

