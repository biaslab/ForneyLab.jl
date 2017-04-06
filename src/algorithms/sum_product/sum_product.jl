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
        entry.msg_update_rule = SumProductRule{typeof(entry.interface.node)}
    end

    inferUpdateRules!(schedule)

    return schedule
end

function inferUpdateRule!(  entry::ScheduleEntry,
                            ::Type{SumProductRule},
                            inferred_outbound_types::Dict{Interface, DataType})
    # Collect inbound types
    inbound_types = DataType[]
    for node_interface in entry.node.interfaces
        if is(node_interface, entry.interface)
            push!(inbound_types, Void)
        else
            push!(inbound_types, inferred_outbound_types[interface.partner])
        end
    end

    # Find applicable rule(s)
    applicable_rules = DataType[]
    for rule in leaftypes(entry.msg_update_rule)
        if isApplicable(rule, inbound_types)
            push!(applicable_rules, rule)
        end
    end

    # Select and set applicable rule
    if isempty(applicable_rules)
        error("No applicable msg update rule for $(entry) with inbound types $(inbound_types)")
    elseif length(applicable_rules) > 1
        error("Multiple applicable msg update rules for $(entry) with inbound types $(inbound_types)")
    else
        entry.msg_update_rule = first(applicable_rules)
    end

    return entry
end