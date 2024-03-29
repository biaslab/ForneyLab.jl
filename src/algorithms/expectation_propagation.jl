export
ExpectationPropagationRule,
@expectationPropagationRule

"""
A non-specific expectation propagation update
"""
abstract type ExpectationPropagationRule{factor_type} <: MessageUpdateRule end

messagePassingSchedule(variable::Variable) = messagePassingSchedule([variable])

function inferUpdateRule!(entry::ScheduleEntry,
                          rule_type::Type{T},
                          inferred_outbound_types::Dict
                         ) where T<:ExpectationPropagationRule
    # Collect inbound types
    inbound_types = collectInboundTypes(entry, rule_type, inferred_outbound_types)

    # Find outbound id
    outbound_id = something(findfirst(isequal(entry.interface), entry.interface.node.interfaces), 0)

    # Find applicable rule(s)
    applicable_rules = Type[]
    for rule in leaftypes(entry.message_update_rule)
        if isApplicable(rule, inbound_types, outbound_id)
            push!(applicable_rules, rule)
        end
    end

    # Select and set applicable rule
    entry.message_update_rule = selectApplicableRule(rule_type, entry, inbound_types, applicable_rules)

    return entry
end

function collectInboundTypes(entry::ScheduleEntry,
                             ::Type{T},
                             inferred_outbound_types::Dict
                            ) where T<:ExpectationPropagationRule
    inbound_message_types = Type[]
    for node_interface in entry.interface.node.interfaces
        if isDeltaConstrained(node_interface.partner)
            push!(inbound_message_types, Message{PointMass})
        else
            push!(inbound_message_types, inferred_outbound_types[node_interface.partner])
        end
    end

    return inbound_message_types
end

"""
`@expectationPropagationRule` registers a expectation propagation update
rule by defining the rule type and the corresponding methods for the `outboundType`
and `isApplicable` functions. If no name (type) for the new rule is passed, a
unique name (type) will be generated. Returns the rule type.
"""
macro expectationPropagationRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    outbound_id = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.args[1] == :(=>)) || error("Invalid call to @expectationPropagationRule")

        if arg.args[2].value == :node_type
            node_type = arg.args[3]
        elseif arg.args[2].value == :outbound_type
            outbound_type = arg.args[3]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for ExpectationPropagationRule should be a Message")
        elseif arg.args[2].value == :inbound_types
            inbound_types = arg.args[3]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[2].value == :outbound_id
            outbound_id = arg.args[3]
        elseif arg.args[2].value == :name
            name = arg.args[3]
        else
            error("Unrecognized field $(arg.args[2].value) in call to @expectationPropagationRule")
        end
    end

    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("EP$(node_type)$(msg_types_hash)")
    end

    # Build validators for isApplicable
    input_type_validators = Expr[]

    push!(input_type_validators, :(length(input_types) == $(length(inbound_types.args))))
    for (i, i_type) in enumerate(inbound_types.args)
        if i_type != :Nothing
            # Only validate inbounds required for message update
            push!(input_type_validators, :(ForneyLab.matches(input_types[$i], $i_type)))
        end
    end

    expr = quote
        struct $name <: ExpectationPropagationRule{$node_type} end
        ForneyLab.outboundType(::Type{$name}) = $outbound_type
        ForneyLab.isApplicable(::Type{$name}, input_types::Vector{<:Type}, outbound_id::Int64) = begin
            $(reduce((current, item) -> :($current && $item), input_type_validators, init = :(outbound_id === $outbound_id)))
        end
    end

    return esc(expr)
end

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInbounds(entry::ScheduleEntry, ::Type{T}) where T<:ExpectationPropagationRule
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    inbounds = Any[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if isClamped(inbound_interface)
            # Hard-code outbound message of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, Message))
        else
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        end
    end

    return inbounds
end
