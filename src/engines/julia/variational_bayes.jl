export variationalAlgorithm, freeEnergyAlgorithm

"""
Create a variational algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalAlgorithm(q_factors::Vector{RecognitionFactor}; file::String="", name::String="")
    q_schedule = variationalSchedule(q_factors)
    marginal_schedule = marginalSchedule(q_factors, q_schedule)

    algo = messagePassingAlgorithm(q_schedule, marginal_schedule, file=file, name=name)

    return algo
end
variationalAlgorithm(q_factor::RecognitionFactor; file::String="", name::String="") = variationalAlgorithm([q_factor]; file=file, name=name)
function variationalAlgorithm(q::RecognitionFactorization=currentRecognitionFactorization())
    algos = "begin\n\n"
    for (id, q_factor) in q.recognition_factors
        algos *= variationalAlgorithm(q_factor, name="$(id)")
        algos *= "\n\n"
    end
    algos *= "\nend # block"
    
    return algos
end

"""
Construct argument code for naive VB updates
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) where T<:NaiveVariationalRule = collectNaiveVariationalNodeInbounds(entry.interface.node, entry, interface_to_msg_idx)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectNaiveVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(inbound_interface.node))
        else
            # Collect marginal from marginal dictionary
            push!(inbounds, "marginals[:$(node_interface.edge.variable.id)]")
        end
    end

    return inbounds
end


"""
Construct argument code for structured VB updates
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) where T<:StructuredVariationalRule = collectStructuredVariationalNodeInbounds(entry.interface.node, entry, interface_to_msg_idx)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectStructuredVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge)
    local_cluster_ids = localRecognitionFactorization(entry.interface.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(inbound_interface.node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbounds
end

"""
The `freeEnergyAlgorithm` function accepts a `RecognitionFactorization` and returns
(if possible) Julia code for computing the variational free energy with respect to 
the argument recognition factorization and corresponding `FactorGraph` (model).
"""
function freeEnergyAlgorithm(q=currentRecognitionFactorization(); name::String="")
    # Write evaluation function for free energy
    energy_block = ""
    entropy_block = ""

    for rf in values(q.recognition_factors)
        hasCollider(rf) && error("Cannot construct localized free energy algorithm. Recognition distribution for factor with id :$(rf.id) does not factor according to local graph structure. This is likely due to a conditional dependence in the posterior distribution (see Bishop p.485). Consider wrapping conditionally dependent variables in a composite node.")
    end

    for node in sort(collect(values(q.graph.nodes)))
        if !isa(node, DeltaFactor) # Non-deterministic factor, add to free energy functional
            # Construct average energy term
            node_str = replace(string(typeof(node)), "ForneyLab." => "") # Remove module prefixes
            inbounds = collectAverageEnergyInbounds(node)
            inbounds_str = join(inbounds, ", ")
            energy_block *= "F += averageEnergy($node_str, $inbounds_str)\n"

            # Construct differential entropy term
            outbound_interface = node.interfaces[1]
            outbound_partner = ultimatePartner(outbound_interface)
            if !(outbound_partner == nothing) && !isa(outbound_partner.node, Clamp) # Differential entropy is required
                dict = q.node_edge_to_cluster
                if haskey(dict, (node, outbound_interface.edge)) # Outbound edge is part of a cluster
                    inbounds = collectConditionalDifferentialEntropyInbounds(node) # Collect conditioning terms for conditional differential entropy
                    inbounds_str = join(inbounds, ", ")
                    entropy_block *= "F -= conditionalDifferentialEntropy($inbounds_str)\n"
                else
                    marginal_idx = outbound_interface.edge.variable.id
                    entropy_block *= "F -= differentialEntropy(marginals[:$marginal_idx])\n"
                end
            end        
        end
    end

    # Combine blocks
    code = "function freeEnergy$(name)(data::Dict, marginals::Dict)\n\n"
    code *= "F = 0.0\n\n"
    code *= energy_block*"\n"*entropy_block
    code *= "\nreturn F\n\n" 
    code *= "end"

    return code
end

function collectAverageEnergyInbounds(node::FactorNode)
    inbounds = String[]

    local_cluster_ids = localRecognitionFactorization(node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(inbound_interface.node))
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbounds
end

function collectConditionalDifferentialEntropyInbounds(node::FactorNode)
    inbounds = String[]

    outbound_edge = node.interfaces[1].edge
    dict = current_recognition_factorization.node_edge_to_cluster
    cluster = dict[(node, outbound_edge)]

    push!(inbounds, "marginals[:$(cluster.id)]") # Add joint term to inbounds

    # Add conditioning terms to inbounds
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)

        if !(node_interface.edge in cluster.edges)
            # Only collect conditioning variables that are part of the cluster
            continue
        elseif (node_interface.edge == outbound_edge)
            # Skip the outbound edge, whose variable is not part of the conditioning term
            continue
        elseif (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in inbounds
            push!(inbounds, marginalString(inbound_interface.node))
        else
            marginal_idx = node_interface.edge.variable.id
            push!(inbounds, "marginals[:$marginal_idx]")
        end
    end

    return inbounds
end