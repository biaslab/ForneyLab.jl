export expectationPropagationAlgorithm, variationalExpectationPropagationAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function expectationPropagationAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = expectationPropagationSchedule(variables)
    marginal_schedule = marginalSchedule(variables)
    
    algo = "begin\n\n"
    algo *= messagePassingAlgorithm(schedule, marginal_schedule, file=file, name=name)
    algo *= "\n\nend # block"

    return algo
end
expectationPropagationAlgorithm(variable::Variable; file::String="", name::String="") = expectationPropagationAlgorithm([variable], file=file, name=name)

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(q_factor::RecognitionFactor; file::String="", name::String="")
    q_schedule = variationalExpectationPropagationSchedule(q_factor)
    marginal_schedule = marginalSchedule(q_factor, q_schedule)

    algo = messagePassingAlgorithm(q_schedule, marginal_schedule, file=file, name=name)

    return algo
end
function variationalExpectationPropagationAlgorithm(q::RecognitionFactorization=currentRecognitionFactorization())
    algos = "begin\n\n"
    for (id, q_factor) in q.recognition_factors
        algos *= variationalExpectationPropagationAlgorithm(q_factor, name="$(id)")
        algos *= "\n\n"
    end
    algos *= "\nend # block"
    
    return algos
end

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) where T<:ExpectationPropagationRule
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
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