export
NaiveVariationalRule,
@naiveVariationalRule

"""
A non-specific naive variational update
"""
abstract type NaiveVariationalRule{factor_type} <: MessageUpdateRule end

"""
Infer the update rule that computes the message for `entry`, as dependent on the inbound types
"""
function inferUpdateRule!(  entry::ScheduleEntry,
                            rule_type::Type{T},
                            inferred_outbound_types::Dict{Interface, Type}) where T<:NaiveVariationalRule
    # Collect inbound types
    inbound_types = collectInboundTypes(entry, rule_type, inferred_outbound_types)

    # Find applicable rule(s)
    applicable_rules = Type[]
    for rule in leaftypes(entry.message_update_rule)
        if isApplicable(rule, inbound_types)
            push!(applicable_rules, rule)
        end
    end

    # Select and set applicable rule
    if isempty(applicable_rules)
        error("No applicable $(rule_type) update for $(typeof(entry.interface.node)) node with inbound types: $(join(inbound_types, ", "))")
    elseif length(applicable_rules) > 1
        error("Multiple applicable $(rule_type) updates for $(typeof(entry.interface.node)) node with inbound types: $(join(inbound_types, ", "))")
    else
        entry.message_update_rule = first(applicable_rules)
    end

    return entry
end

"""
Find the inbound types that are required to compute the message for `entry`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(   entry::ScheduleEntry,
                                ::Type{T},
                                inferred_outbound_types::Dict{Interface, Type}) where T<:NaiveVariationalRule
    inbound_types = Type[]
    for node_interface in entry.interface.node.interfaces
        if node_interface === entry.interface
            push!(inbound_types, Nothing)
        else
            # Edge is external, accept marginal
            push!(inbound_types, ProbabilityDistribution)
        end
    end

    return inbound_types
end

"""
`@naiveVariationalRule` registers a variational update rule for the naive
(mean-field) factorization by defining the rule type and the corresponding
methods for the `outboundType` and `isApplicable` functions. If no name (type)
for the new rule is passed, a unique name (type) will be generated. Returns the
rule type.
"""
macro naiveVariationalRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.args[1] == :(=>)) || error("Invalid call to @naiveVariationalRule")

        if arg.args[2].value == :node_type
            node_type = arg.args[3]
        elseif arg.args[2].value == :outbound_type
            outbound_type = arg.args[3]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for VariationalRule should be a Message")
        elseif arg.args[2].value == :inbound_types
            inbound_types = arg.args[3]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[2].value == :name
            name = arg.args[3]
        else
            error("Unrecognized field $(arg.args[2].value) in call to @naiveVariationalRule")
        end
    end

    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("VB$(node_type)$(msg_types_hash)")
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
        struct $name <: NaiveVariationalRule{$node_type} end
        ForneyLab.outboundType(::Type{$name}) = $outbound_type
        ForneyLab.isApplicable(::Type{$name}, input_types::Vector{<:Type}) = begin
            $(reduce((current, item) -> :($current && $item), input_type_validators, init = :true))
        end
    end

    return esc(expr)
end

"""
Construct argument code for naive VB updates
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}) where T<:NaiveVariationalRule = collectNaiveVariationalNodeInbounds(entry.interface.node, entry)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectNaiveVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry)
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry

    inbounds = Any[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface === entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, nothing)
        elseif isClamped(inbound_interface)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, ProbabilityDistribution))
        else
            # Collect entry from marginal schedule
            push!(inbounds, target_to_marginal_entry[node_interface.edge.variable])
        end
    end

    return inbounds
end
