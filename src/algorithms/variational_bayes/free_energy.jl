"""
The `assembleFreeEnergy` function accepts an `Algorithm` and populates
required fields for computing the variational free energy.
"""
function assembleFreeEnergy!(algo=currentAlgorithm())
    # Find counting numbers for energies and entropies
    assembleCountingNumbers!(algo)

    # Convert energy counting numbers to energy inbounds
    average_energies_vect = Vector{Dict{Symbol, Any}}()
    entropies_vect = Vector{Dict{Symbol, Any}}()

    for (node, cnt) in algo.energy_counting_numbers
        if cnt != 0
            average_energy = Dict{Symbol, Any}(:counting_number => cnt,
                                               :node => typeof(node),
                                               :inbounds => collectAverageEnergyInbounds(node))
            push!(average_energies_vect, average_energy)
        end
    end

    # Convert entropy counting numbers to entropy inbounds
    for (target, cnt) in algo.entropy_counting_numbers
        if cnt != 0
            entropy = Dict{Symbol, Any}(:counting_number => cnt,
                                        :inbound => algo.target_to_marginal_entry[target])
            push!(entropies_vect, entropy)
        end
    end

    algo.average_energies = average_energies_vect
    algo.entropies = entropies_vect
    
    return algo    
end

"""
The `assembleCountingNumbers` function accepts an `Algorithm` and
populates the counting numbers for the average energies and entropies.
"""
function assembleCountingNumbers!(algo=currentAlgorithm())
    energy_counting_numbers = Dict{FactorNode, Int64}()
    entropy_counting_numbers = Dict{Union{Variable, Cluster}, Int64}()

    # Collect regions
    internal_edges = Set{Edge}()
    for (id, rf) in algo.recognition_factors
        union!(internal_edges, rf.internal_edges)
    end
    nodes_connected_to_internal_edges = nodes(internal_edges)
    
    # Iterate over large regions
    for node in nodes_connected_to_internal_edges
        if !isa(node, DeltaFactor) # Node is stochastic
            increase!(energy_counting_numbers, node, 1) # Count average energy
            for target in unique!(localClusters(node)) # Collect all unique clusters/variables around node
                if first(target.edges) in internal_edges # Cluster/variable is internal to a recognition factor
                    increase!(entropy_counting_numbers, target, 1) # Count (joint) entropy
                end
            end
        elseif isa(node, Equality)
            increase!(entropy_counting_numbers, node.i[1].edge.variable, 1) # Count univariate entropy
        elseif !isa(node, Clamp) # Node is deterministic and not clamp or equality
            target = cluster(node, node.interfaces[2].edge) # Find cluster/variable of inbound edges
            increase!(entropy_counting_numbers, target, 1) # Count (joint) entropy
        end
    end

    # Iterate over small regions
    for edge in internal_edges
        increase!(entropy_counting_numbers, edge.variable, -1) # Discount univariate entropy
    end

    algo.energy_counting_numbers = energy_counting_numbers
    algo.entropy_counting_numbers = entropy_counting_numbers

    return algo
end

"""
Increase (or decrease) a counting number
"""
function increase!(dict::Dict, key::Any, increase::Number)
    if haskey(dict, key)
        dict[key] += increase
    else
        dict[key] = increase
    end

    return dict
end

function collectAverageEnergyInbounds(node::FactorNode)
    inbounds = Any[]
    local_clusters = localRecognitionFactorization(node)

    recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        node_interface_recognition_factor = recognitionFactor(node_interface.edge)

        if (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif !(node_interface_recognition_factor in recognition_factors)
            # Collect marginal entry from marginal dictionary (if marginal entry is not already accepted)
            target = local_clusters[node_interface_recognition_factor]
            push!(inbounds, current_algorithm.target_to_marginal_entry[target])
        end

        push!(recognition_factors, node_interface_recognition_factor)
    end

    return inbounds
end