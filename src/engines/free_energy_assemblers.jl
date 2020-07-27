"""
The `assembleFreeEnergy` function accepts an `Algorithm` and populates
required fields for computing the variational free energy.
"""
function assembleFreeEnergy!(algo=currentInferenceAlgorithm())
    # Find counting numbers for energies and entropies
    pfz = algo.posterior_factorization
    assembleCountingNumbers!(pfz)

    # Convert energy counting numbers to energy inbounds
    average_energies_vect = Vector{Dict{Symbol, Any}}()
    entropies_vect = Vector{Dict{Symbol, Any}}()

    energy_counting_numbers_vect = [(node, cnt) for (node, cnt) in pfz.energy_counting_numbers] # Order the Dict in a Vector
    for (node, cnt) in sort(energy_counting_numbers_vect)
        if cnt != 0
            average_energy = Dict{Symbol, Any}(:counting_number => cnt,
                                               :node => typeof(node),
                                               :inbounds => collectAverageEnergyInbounds(node))
            push!(average_energies_vect, average_energy)
        end
    end

    # Convert entropy counting numbers to entropy inbounds
    entropy_counting_numbers_vect = [(target, cnt) for (target, cnt) in pfz.entropy_counting_numbers] # Order the Dict in a Vector
    for (target, cnt) in sort(entropy_counting_numbers_vect)
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
function assembleCountingNumbers!(pfz=currentPosteriorFactorization())
    pfz.free_energy_flag || error("Required quantities for free energy evaluation are not computed by the algorithm. Make sure to flag free_energy=true upon algorithm construction to schedule computation of required quantities.")

    energy_counting_numbers = Dict{FactorNode, Int64}()
    entropy_counting_numbers = Dict{Region, Int64}()

    # Collect regions
    internal_edges = Set{Edge}()
    for (_, pf) in pfz.posterior_factors
        union!(internal_edges, pf.internal_edges)
    end
    nodes_connected_to_internal_edges = nodes(internal_edges)
    
    # Iterate over large regions
    for node in nodes_connected_to_internal_edges
        if isa(node, Clamp)
            continue
        elseif !isa(node, DeltaFactor) # Node is stochastic
            increase!(energy_counting_numbers, node, 1) # Count average energy
            for target in unique!(localStochasticRegions(node)) # Collect all unique regions around node
                increase!(entropy_counting_numbers, target, 1) # Count (joint) entropy
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
        increase!(entropy_counting_numbers, edge.variable, -(degree(edge) - 1)) # Discount univariate entropy
    end

    pfz.energy_counting_numbers = energy_counting_numbers
    pfz.entropy_counting_numbers = entropy_counting_numbers

    return pfz
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