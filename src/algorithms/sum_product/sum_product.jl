export
SumProductRule,
sumProductSchedule,
@sumProductRule

abstract type SumProductRule{factor_type} <: MessageUpdateRule end

"""
sumProductSchedule() generates a sum-product message passing schedule that computes the
marginals for each of the argument variables.
"""
function sumProductSchedule(variables::Vector{Variable})
    # Generate a feasible summary propagation schedule
    schedule = summaryPropagationSchedule(variables)

    # Assign the sum-product update rule to each of the schedule entries
    for entry in schedule
        entry.msg_update_rule = SumProductRule{typeof(entry.interface.node)}
    end

    inferUpdateRules!(schedule)

    return schedule
end

sumProductSchedule(variable::Variable) = sumProductSchedule([variable])

function inferUpdateRule!{T<:SumProductRule}(   entry::ScheduleEntry,
                                                rule_type::Type{T},
                                                inferred_outbound_types::Dict{Interface, DataType})
    # Collect inbound types
    inbound_types = collectInboundTypes(entry, rule_type, inferred_outbound_types)

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
        error("Multiple applicable msg update rules for $(entry) with inbound types $(inbound_types): $(applicable_rules)")
    else
        entry.msg_update_rule = first(applicable_rules)
    end

    return entry
end

function collectInboundTypes{T<:SumProductRule}(entry::ScheduleEntry,
                                                ::Type{T},
                                                inferred_outbound_types::Dict{Interface, DataType})
    inbound_message_types = DataType[]
    for node_interface in entry.interface.node.interfaces
        if node_interface == entry.interface
            push!(inbound_message_types, Void)
        elseif (node_interface.partner != nothing) && isa(node_interface.partner.node, Clamp)
            push!(inbound_message_types, Message{PointMass})
        else
            push!(inbound_message_types, inferred_outbound_types[node_interface.partner])
        end
    end

    return inbound_message_types
end

"""
@sumProductRule registers a sum-product update rule by defining the rule type
and the corresponding methods for the outboundType and isApplicable functions.
If no name (type) for the new rule is passed, a unique name (type) will be
generated. Returns the rule type.
"""
macro sumProductRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
     
        (arg.args[1] == :(=>)) || error("Invalid call to @sumProductRule")

        if arg.args[2].args[1] == :node_type
            node_type = arg.args[3]
        elseif arg.args[2].args[1] == :outbound_type
            outbound_type = arg.args[3]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for SumProductRule should be a Message")
        elseif arg.args[2].args[1] == :inbound_types
            inbound_types = arg.args[3]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[2].args[1] == :name
            name = arg.args[3]
        else
            error("Unrecognized field $(arg.args[2].args[1]) in call to @sumProductRule")
        end
    end
    
    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("SP$(node_type)$(msg_types_hash)")
    end

    # Build validators for isApplicable
    input_type_validators = String[]
    for (i, i_type) in enumerate(inbound_types.args)
        if i_type != :Void
            # Only validate inbounds required for message update
            push!(input_type_validators, "ForneyLab.matches(input_types[$i], $i_type)")
        end
    end

    expr = parse("""
        begin
            type $name <: SumProductRule{$node_type} end
            ForneyLab.outboundType(::Type{$name}) = $outbound_type
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{DataType}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end