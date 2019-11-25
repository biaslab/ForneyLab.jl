export expectationPropagationAlgorithm, variationalExpectationPropagationAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function expectationPropagationAlgorithm(variables::Vector{Variable}; name::String="")
    schedule = expectationPropagationSchedule(variables)
    marginal_schedule = marginalSchedule(variables)
    
    # Build (empty) recognition factor datastructure
    rf_dict = messagePassingAlgorithm(schedule, marginal_schedule)
    rf_dict[:id] = Symbol("")

    # Build algorithm datastructure
    algo_dict = Dict{Symbol, Any}(:name => name,
                                  :recognition_factors => [rf_dict])

    algo_str = algorithmString(algo_dict)
    
    return algo_str
end
expectationPropagationAlgorithm(variable::Variable; name::String="") = expectationPropagationAlgorithm([variable], name=name)

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(q::RecognitionFactorization=currentRecognitionFactorization())
    recognition_factors_vect = Vector{Dict{Symbol, Any}}(undef, length(q.recognition_factors))
    algo_dict = Dict{Symbol, Any}(:name => name,
                                  :recognition_factors => recognition_factors_vect)

    for (i, (id, q_factor)) in enumerate(q.recognition_factors)
        schedule = variationalExpectationPropagationSchedule(q_factor)
        marginal_schedule = marginalSchedule(q_factor, schedule)

        # Populate algorithm datastructure
        algo_dict[:recognition_factors][i] = messagePassingAlgorithm(schedule, marginal_schedule)
        algo_dict[:recognition_factors][i][:id] = q_factor.id
    end

    algo_str = algorithmString(algo_dict)

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
            push!(inbounds, messageDict(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, Dict{Symbol, Any}(:schedule_index => inbound_idx))
        end
    end

    return inbounds
end