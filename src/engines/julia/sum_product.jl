export sumProductAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function sumProductAlgorithm(variables::Vector{Variable}, rfz::RecognitionFactorization=currentRecognitionFactorization())
    # Initialize a container recognition factor
    rf = RecognitionFactor(rfz, id=Symbol(""))
    schedule = sumProductSchedule(variables)
    rf.schedule = condense(flatten(schedule)) # Inline all internal message passing and remove clamp node entries
    rf.marginal_table = marginalTable(variables)

    assembleAlgorithm!(rf)
    algo_str = algorithmString(rfz)
    
    return algo_str
end
sumProductAlgorithm(variable::Variable, rfz::RecognitionFactorization=currentRecognitionFactorization()) = sumProductAlgorithm([variable], rfz)

"""
Collect and construct SP update code for each inbound.
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}, interface_to_schedule_entry::Dict, ::Dict) where T<:SumProductRule = collectSumProductNodeInbounds(entry.interface.node, entry, interface_to_schedule_entry)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectSumProductNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_schedule_entry::Dict)
    inbounds = Any[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface === entry.interface
            # Ignore inbound message on outbound interface
            push!(inbounds, nothing)
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    return inbounds
end
