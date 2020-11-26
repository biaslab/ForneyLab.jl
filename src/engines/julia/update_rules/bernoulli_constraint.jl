export ruleSPBernoulliConstraintOut


#-------------
# Update rules
#-------------

function ruleSPBernoulliConstraintOut(msg_out::Message{Bernoulli, Univariate}, p::Float64)
    p_vec = [p, 1-p]
    p_out_vec = [msg_out.dist.params[:p], 1-msg_out.dist.params[:p]]
    p_div = clamp.(p_vec ./ p_out_vec, tiny, huge) # Soften vanishing probabilities
    norm = sum(p_div)
    p_norm_vec = p_div./norm

    return Message(Univariate, Bernoulli, p=p_norm_vec[1])
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::BernoulliConstraint, entry::ScheduleEntry)
    inbounds = Any[]

    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    # Push function and value to calling signature; these need to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:p => node.p,
                                      :keyword => false))
    return inbounds
end