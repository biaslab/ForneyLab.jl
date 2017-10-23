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
        internal_edges = extend(edges(variables))
        self = new(id, variables, internal_edges)
        rfz.recognition_factors[id] = self # Register new factor with recognition factorization

        return self 
    end
end

RecognitionFactor(variable::Variable; id=generateId(RecognitionFactor)) = RecognitionFactor(Set([variable]), id=id)
RecognitionFactor(variables::Vector{Variable}; id=generateId(RecognitionFactor)) = RecognitionFactor(Set(variables), id=id)

"""
factor() returns a (list of) recognition factor(s) that can be used for convenient scheduling
"""
function factor(variables::Set{Variable}; factorization=:structured, id=generateId(RecognitionFactor))
    if factorization == :naive
        # Factorize the subgraph in deterministically-connected clusters.
        # This is convenient when constructing schedules, because
        # several clusters can be grouped under one schedule/algorithm.

        cluster_sig_to_variables_dict = Dict{UInt64, Vector{Variable}}()
        # Find which variables belong to the same cluster
        for variable in sort(collect(variables))
            cluster_sig = hash(extend(edges(variable))) # Generate a unique cluster signature
            if haskey(cluster_sig_to_variables_dict, cluster_sig)
                # Cluster is already registered, push variable to registered cluster
                push!(cluster_sig_to_variables_dict[cluster_sig], variable)
                # TODO: this might lead to problems when computing the free energy
            else
                # Register new cluster
                cluster_sig_to_variables_dict[cluster_sig] = [variable]
            end
        end

        # For each unique cluster, construct a recognition factor
        recognition_factors = RecognitionFactor[]
        for (i, cluster_sig) in enumerate(sort(collect(keys(cluster_sig_to_variables_dict))))
            rf = RecognitionFactor(cluster_sig_to_variables_dict[cluster_sig], id=id*:_*i)
            push!(recognition_factors, rf)
        end
        return recognition_factors
    elseif factorization == :structured
        # 
        return RecognitionFactor(variables, id=id)
    else
        error("Unknown factorization keyword $(factorization), choose `:naive` or `:structured` instead")
    end
end
factor(variable::Variable; factorization=:structured, id=generateId(RecognitionFactor)) = factor(Set([variable]), factorization=factorization, id=id)
factor(variables::Vector{Variable}; factorization=:structured, id=generateId(RecognitionFactor)) = factor(Set(variables), factorization=factorization, id=id)

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