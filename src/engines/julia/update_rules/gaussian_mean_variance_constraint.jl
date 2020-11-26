export ruleSPGaussianMeanVarianceConstraintOut


#-------------
# Update rules
#-------------

function ruleSPGaussianMeanVarianceConstraintOut(msg_out::Message{<:Gaussian, V}, 
                                                 m::Union{Float64, AbstractVector}, 
                                                 v::Union{Float64, AbstractMatrix}) where V<:VariateType

    (xi_out, w_out) = unsafeWeightedMeanPrecision(msg_out.dist)
    w_q = cholinv(v)
    xi_q = w_q*m

    return Message(V, GaussianWeightedMeanPrecision, xi=xi_q-xi_out, w=w_q-w_out)
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