export ruleVBCVIOutVD, ruleVBCVIIn1MV, ruleVBCVIOutVDX, ruleVBCVIInX

function ruleVBCVIIn1MV(node_id::Symbol,
                        msg_out::ProbabilityDistribution,
                        msg_in::Message{<:FactorNode, <:VariateType})

    @show msg_in
    @show msg_out
    msg_in
end

function ruleVBCVIOutVD(node_id::Symbol,
                        msg_out::Any,
                        msg_in::ProbabilityDistribution)

    thenode = currentGraph().nodes[node_id]

    samples = thenode.g.(sample(msg_in, thenode.num_samples))
    weights = ones(thenode.num_samples)/thenode.num_samples

    if length(samples[1]) == 1
        variate = Univariate
    else
        variate = Multivariate
    end

    q=ProbabilityDistribution(variate, SampleList, s=samples, w=weights)
    q.params[:entropy] = 0

    return Message(variate,SetSampleList,q=q,node_id=node_id)

end

function ruleVBCVIInX(node_id::Symbol,
                      inx::Int64,
                      msg_out::ProbabilityDistribution,
                      msgs_in::Vararg{Union{Message,ProbabilityDistribution}})

    @show msgs_in
    msgs_in[inx]
end

function ruleVBCVIOutVDX(node_id::Symbol,
                         msg_out::Any,
                         msgs_in::Vararg{ProbabilityDistribution})

    thenode = currentGraph().nodes[node_id]

    samples_in = [sample(msg_in, thenode.num_samples) for msg_in in msgs_in]

    samples = thenode.g.(samples_in...)
    weights = ones(thenode.num_samples)/thenode.num_samples

    if length(samples[1]) == 1
        variate = Univariate
    else
        variate = Multivariate
    end

    q=ProbabilityDistribution(variate, SampleList, s=samples, w=weights)
    q.params[:entropy] = 0

    return Message(variate,SetSampleList,q=q,node_id=node_id)

end


#---------------------------
# Custom inbounds collectors
#---------------------------

# function collectStructuredVariationalNodeInbounds(node::CVI, entry::ScheduleEntry)
#     interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
#     target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry
#
#     inbounds = Any[]
#     entry_posterior_factor = posteriorFactor(entry.interface.edge)
#     local_edge_to_region = localEdgeToRegion(entry.interface.node)
#
#     push!(inbounds, node.id)
#
#     multi_in = (length(node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
#     inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound
#
#     if (inx > 0) && multi_in # Multi-inbound backward rule
#         push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier
#                                           :keyword => false))
#     end
#
#     encountered_posterior_factors = Union{PosteriorFactor, Edge}[] # Keep track of encountered posterior factors
#     for node_interface in node.interfaces
#         inbound_interface = ultimatePartner(node_interface)
#         current_posterior_factor = posteriorFactor(node_interface.edge)
#
#         if (node_interface == entry.interface != node.interfaces[1])
#             # Collect the incoming message
#             push!(inbounds, interface_to_schedule_entry[inbound_interface])
#         elseif node_interface === entry.interface
#             # Ignore marginal of outbound edge
#             push!(inbounds, nothing)
#         elseif isClamped(inbound_interface)
#             # Hard-code marginal of constant node in schedule
#             push!(inbounds, assembleClamp!(inbound_interface.node, ProbabilityDistribution))
#         elseif current_posterior_factor === entry_posterior_factor
#             # Collect message from previous result
#             push!(inbounds, interface_to_schedule_entry[inbound_interface])
#         elseif !(current_posterior_factor in encountered_posterior_factors)
#             # Collect marginal from marginal dictionary (if marginal is not already accepted)
#             target = local_edge_to_region[node_interface.edge]
#             push!(inbounds, target_to_marginal_entry[target])
#         end
#
#         push!(encountered_posterior_factors, current_posterior_factor)
#     end
#
#     return inbounds
# end

function collectNaiveVariationalNodeInbounds(node::CVI, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry

    inbounds = Any[]

    push!(inbounds, node.id)

    multi_in = (length(node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
    inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound

    if (inx > 0) && multi_in # Multi-inbound backward rule
        push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier
                                          :keyword => false))
    end

    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)

        if (node_interface == entry.interface != node.interfaces[1])
            # Collect the incoming message
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif node_interface === entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, nothing)
        elseif isClamped(inbound_interface)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, ProbabilityDistribution))
        else
            # Collect entry from marginal schedule
            push!(inbounds, target_to_marginal_entry[node_interface.edge.variable])
        end
    end

    return inbounds
end
