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

"""
Construct argument code for naive VB updates
"""
collectInbounds{T<:NaiveVariationalRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) = collectNaiveVariationalNodeInbounds(entry.interface.node, entry, interface_to_msg_idx)

"""
Construct the inbound code that computes the message for `entry`.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectInbounds{T<:StructuredVariationalRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge)
    local_cluster_ids = localRecognitionFactorization(entry.interface.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.interface.node.interfaces
        inbound_interface = node_interface.partner
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
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
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectNaiveVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = node_interface.partner
        partner_node = inbound_interface.node
        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        else
            # Collect marginal from marginal dictionary
            push!(inbounds, "marginals[:$(node_interface.edge.variable.id)]")
        end
    end

    return inbounds
end

"""
Collect marginals associated with all edges connected to `node`
"""
# TODO: use collectInbounds for marginals?
# function collectMarginals(node::FactorNode)
#     # Collect marginals
#     inbound_marginals = String[]
#     for interface in node.interfaces
#         partner_node = interface.partner.node
#         if isa(partner_node, Clamp)
#             # Hard-code marginal of constant node
#             push!(inbound_marginals, marginalString(partner_node))
#         else
#             # Collect marginal from marginal dictionary
#             push!(inbound_marginals, "marginals[:$(interface.edge.variable.id)]")
#         end
#     end

#     return inbound_marginals
# end

# TODO: fix for structured factorization
# function freeEnergyAlgorithm(recognition_factors::Vector{RecognitionFactor}=collect(values(current_recognition_factorization.recognition_factors)))
#     # Collect nodes connected to external edges
#     nodes_connected_to_external_edges = Set{FactorNode}()
#     for rf in recognition_factors
#         union!(nodes_connected_to_external_edges, nodesConnectedToExternalEdges(rf))
#     end

#     # Write evaluation function for free energy
#     energy_block = ""
#     entropy_block = ""
#     for node in sort(collect(nodes_connected_to_external_edges))
#         # Average energy
#         inbounds = collectMarginals(node)
#         inbounds_str = join(inbounds, ", ")
#         node_str = replace(string(typeof(node)),"ForneyLab.", "") # Remove module prefixes
#         energy_block *= "F += averageEnergy($(node_str), $(inbounds_str))\n"
        
#         # Differential entropy
#         if !(node.interfaces[1].partner == nothing) && !isa(node.interfaces[1].partner.node, Clamp)
#             entropy_block *= "F -= differentialEntropy(marginals[:$(node.interfaces[1].edge.variable.id)])\n"
#         end
#     end

#     # Combine blocks
#     code = "function freeEnergy(data::Dict, marginals::Dict)\n\n"
#     code *= "F = 0.0\n\n"
#     code *= energy_block*"\n"entropy_block
#     code *= "\nreturn F\n" 
#     code *= "end"

#     return code
# end