export ruleSPCategoricalConstraintOut


#-------------
# Update rules
#-------------

function ruleSPCategoricalConstraintOut(msg_out::Message{Categorical, Univariate}, p::Vector{Float64})
    p_div = clamp.(p ./ msg_out.dist.params[:p], tiny, huge) # Soften vanishing probabilities
    norm = sum(p_div)

    return Message(Univariate, Categorical, p=p_div./norm)
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::CategoricalConstraint, entry::ScheduleEntry)
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