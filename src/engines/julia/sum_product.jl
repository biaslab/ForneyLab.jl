export sumProductAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function sumProductAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = sumProductSchedule(variables)
    algo = messagePassingAlgorithm(schedule, variables, file=file, name=name)

    return algo
end
sumProductAlgorithm(variable::Variable; file::String="", name::String="") = sumProductAlgorithm([variable], file=file, name=name)

"""
Collect and construct SP update code for each inbound.
"""
function collectInbounds{T<:SumProductRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = node_interface.partner
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