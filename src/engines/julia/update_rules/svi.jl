export ruleSPSVIIn1MN, ruleSPSVIOutNM

#----------------------------------
# SPSVIIn1MN must be called before SPSVIOutNM in the message passing schedular.
#----------------------------------

function ruleSPSVIIn1MN(node_id::Symbol,
                        msg_out::Message{<:FactorFunction, <:VariateType},
                        msg_in::Message{<:FactorFunction, <:VariateType})

    thenode = currentGraph().nodes[node_id]
    mod(thenode.dataset_size/thenode.batch_size,1) == 0 || error("Data size must be integer multiple of the batch size.") # must be in isApplicable
    scale = Int(thenode.dataset_size/thenode.batch_size)

    # Calculate the message in such a way that the variational distribution of the global variable
    # will match the stochastic update of variational distribution.
    η = deepcopy(naturalParams(thenode.q))
    q_ = deepcopy(msg_in.dist)
    m_out_dist = deepcopy(msg_out.dist)
    for _=1:scale q_ = prod!(q_,m_out_dist) end
    ∇ = η .- naturalParams(q_) # natural gradient
    update!(thenode.opt,η,∇) # η is updated
    λ = η .- naturalParams(msg_in.dist) # natural parameters for the message

    # Update q stored in the SVI node
    thenode.q = standardDist(thenode.q,η)

    return standardMessage(thenode.q,λ)

end

function ruleSPSVIOutNM(node_id::Symbol,
                        msg_out::Nothing,
                        msg_in::Message{<:FactorFunction, V}) where V<:VariateType

    thenode = currentGraph().nodes[node_id]

    # The out variable (the mirror variable) is not updated yet! It will be used in calculation of other variational distributions.
    # Set the variational distribution of out variable to q that is not updated yet!
    m_out = Message(V,SetProbDist,q=deepcopy(thenode.q),message=msg_in,node_id=node_id)

    return m_out

end

#---------------------------
# Custom inbounds collectors
#---------------------------

function collectSumProductNodeInbounds(node::SVI, entry::ScheduleEntry)
    inbounds = Any[]

    push!(inbounds, node.id)

    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if (node_interface == entry.interface != node.interfaces[1])
            # Collect the message
            #haskey(interface_to_schedule_entry, inbound_interface) || error("The nonlinear node's backward rule uses the incoming message on the input edge to determine the approximation point. Try altering the variable order in the scheduler to first perform a forward pass.")
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbounds, nothing)
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    return inbounds
end
