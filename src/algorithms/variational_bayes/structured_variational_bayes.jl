export
StructuredVariationalRule,
@structuredVariationalRule

abstract type StructuredVariationalRule{factor_type} <: MessageUpdateRule end

"""
Infer the update rule that computes the message for `entry`, as dependent on the inbound types
"""
function inferUpdateRule!(entry::ScheduleEntry,
                          rule_type::Type{T},
                          inferred_outbound_types::Dict{Interface, Type}
                         ) where T<:StructuredVariationalRule
    # Collect inbound types
    inbound_types = collectInboundTypes(entry, rule_type, inferred_outbound_types)
    
    # Find applicable rule(s)
    applicable_rules = Type[]
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

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(entry::ScheduleEntry,
                             ::Type{T},
                             inferred_outbound_types::Dict{Interface, Type}
                            ) where T<:StructuredVariationalRule
    inbound_types = Type[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge) # Recognition factor id for outbound edge
    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.interface.node.interfaces
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if node_interface == entry.interface
            push!(inbound_types, Nothing)
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Edge is internal, accept message
            push!(inbound_types, inferred_outbound_types[node_interface.partner])
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Edge is external, accept marginal (if marginal is not already accepted)
            push!(inbound_types, ProbabilityDistribution) 
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbound_types
end

""" 
`@structuredVariationalRule` registers a variational update rule for the
structured factorization by defining the rule type and the corresponding methods
for the `outboundType` and `isApplicable` functions. If no name (type) for the
new rule is passed, a unique name (type) will be generated. Returns the rule
type.
""" 
macro structuredVariationalRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.args[1] == :(=>)) || error("Invalid call to @structuredVariationalRule")

        if arg.args[2].value == :node_type
            node_type = arg.args[3]
        elseif arg.args[2].value == :outbound_type
            outbound_type = arg.args[3]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for StructuredVariationalRule should be a Message")
        elseif arg.args[2].value == :inbound_types
            inbound_types = arg.args[3]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[2].value == :name
            name = arg.args[3]
        else
            error("Unrecognized field $(arg.args[2].value) in call to @structuredVariationalRule")
        end
    end

    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("SVB$(node_type)$(msg_types_hash)")
    end

    # Build validators for isApplicable
    input_type_validators = String[]
    for (i, i_type) in enumerate(inbound_types.args)
        if i_type != :Nothing
            # Only validate inbounds required for message update
            push!(input_type_validators, "ForneyLab.matches(input_types[$i], $i_type)")
        end
    end

    expr = parse("""
        begin
            mutable struct $name <: StructuredVariationalRule{$node_type} end
            ForneyLab.outboundType(::Type{$name}) = $outbound_type
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{<:Type}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end