export expectationPropagationAlgorithm, variationalExpectationPropagationAlgorithm

"""
Create an expectation propagation algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function expectationPropagationAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = expectationPropagationSchedule(variables)
    marginal_schedule = marginalSchedule(variables)

    algo = messagePassingAlgorithm(schedule, marginal_schedule, file=file, name=name)

    return algo
end
expectationPropagationAlgorithm(variable::Variable; file::String="", name::String="") = expectationPropagationAlgorithm([variable], file=file, name=name)

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(q_factor::RecognitionFactor; file::String="", name::String="")
    q_schedule = variationalExpectationPropagationSchedule(q_factor)
    marginal_schedule = marginalSchedule(q_factor)

    algo = messagePassingAlgorithm(q_schedule, marginal_schedule, file=file, name=name)

    return algo
end

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInbounds{T<:ExpectationPropagationRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = node_interface.partner
        if isa(inbound_interface.node, Clamp)
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