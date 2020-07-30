export ruleSPPointMassConstraintOutG


#-------------
# Update rules
#-------------

ruleSPPointMassConstraintOutG(msg_out::Message{<:Gaussian, Univariate}) = Message(Univariate, PointMass, m=mode(msg_out.dist))


#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::PointMassConstraint, entry::ScheduleEntry)
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

    return inbounds
end