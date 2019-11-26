export expectationPropagationAlgorithm, variationalExpectationPropagationAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function expectationPropagationAlgorithm(variables::Vector{Variable}, rfz::RecognitionFactorization=currentRecognitionFactorization())
    schedule = expectationPropagationSchedule(variables)
    marginal_schedule = marginalSchedule(variables)
    
    # Assemble algorithm in a container recognition factor
    rf = RecognitionFactor(rfz, id=Symbol(""))
    assembleAlgorithm!(rf, schedule, marginal_schedule)
    algo_str = algorithmString(rfz)
    
    return algo_str
end
expectationPropagationAlgorithm(variable::Variable, rfz::RecognitionFactorization=currentRecognitionFactorization()) = expectationPropagationAlgorithm([variable], rfz)

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(rfz::RecognitionFactorization=currentRecognitionFactorization())
    for (id, rf) in rfz.recognition_factors
        schedule = variationalExpectationPropagationSchedule(rf)
        marginal_schedule = marginalSchedule(rf, schedule)
        assembleAlgorithm!(rf, schedule, marginal_schedule)
    end

    algo_str = algorithmString(rfz)

    return algo_str
end

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) where T<:ExpectationPropagationRule
    inbounds = Dict{Symbol, Any}[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleMessageInbound(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, Dict{Symbol, Any}(:schedule_index => inbound_idx))
        end
    end

    return inbounds
end