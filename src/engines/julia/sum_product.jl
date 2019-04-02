export sumProductAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function sumProductAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = sumProductSchedule(variables)
    marginal_schedule = marginalSchedule(variables)

    algo = "begin\n\n"
    algo *= messagePassingAlgorithm(schedule, marginal_schedule, file=file, name=name)
    algo *= "\n\nend # block"
    
    return algo
end
sumProductAlgorithm(variable::Variable; file::String="", name::String="") = sumProductAlgorithm([variable], file=file, name=name)

"""
Collect and construct SP update code for each inbound.
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) where T<:SumProductRule = collectSumProductNodeInbounds(entry.interface.node, entry, interface_to_msg_idx)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectSumProductNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbound_messages, "nothing")
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    return inbound_messages
end
