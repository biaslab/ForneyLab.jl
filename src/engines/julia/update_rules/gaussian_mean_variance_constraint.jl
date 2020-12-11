export ruleSPGaussianMeanVarianceConstraintOut


#-------------
# Update rules
#-------------

function ruleSPGaussianMeanVarianceConstraintOut(msg_out::Message{<:Gaussian, V}, 
                                                 m_f::Union{Float64, AbstractVector}, 
                                                 v_f::Union{Float64, AbstractMatrix}) where V<:VariateType

    (xi_out, w_out) = unsafeWeightedMeanPrecision(msg_out.dist)
    w_f = cholinv(v_f)
    xi_f = w_f*m_f

    return Message(V, GaussianWeightedMeanPrecision, xi=0.5*(xi_f-xi_out), w=0.5*(w_f-w_out))
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::GaussianMeanVarianceConstraint, entry::ScheduleEntry)
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
    push!(inbounds, Dict{Symbol, Any}(:m => node.m,
                                      :keyword => false))
    push!(inbounds, Dict{Symbol, Any}(:v => node.v,
                                      :keyword => false))
    return inbounds
end