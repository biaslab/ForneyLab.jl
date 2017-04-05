export
SumProductRule,
SumProduct

abstract SumProductRule{factor_type} <: MessageUpdateRule

"""
SumProduct() generates a sum-product message passing schedule that computes the
marginals for each of the argument variables.
"""
function SumProduct(variables::Vector{Variable})
    # Generate a feasible summary propagation schedule
    schedule = summaryPropagationSchedule(variables)

    # Assign the sum-product update rule to each of the schedule entries
    for entry in schedule
        entry.msg_update_rule = SumProductRule
    end

    return schedule
end