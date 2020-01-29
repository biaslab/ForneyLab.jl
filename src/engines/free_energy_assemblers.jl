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
    algo.free_energy_flag || error("Required quantities for free energy evaluation are not computed by the algorithm. Make sure to flag free_energy=true upon algorithm construction to schedule computation of required quantities.")

    energy_counting_numbers = Dict{FactorNode, Int64}()
    entropy_counting_numbers = Dict{Region, Int64}()

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
            for target in unique!(localRegions(node)) # Collect all unique regions around node
                if first(target.edges) in internal_edges # Region is internal to a recognition factor
                    increase!(entropy_counting_numbers, target, 1) # Count (joint) entropy
                end
            end
        elseif isa(node, Equality)
            increase!(entropy_counting_numbers, node.i[1].edge.variable, 1) # Count univariate entropy
        elseif length(node.interfaces) >= 2 # Node is deterministic and not equality
            target = region(node, node.interfaces[2].edge) # Find region of inbound edges
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