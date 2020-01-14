export
StructuredVariationalRule,
@structuredVariationalRule

"""
A non-specific structured variational update
"""
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
    for rule in leaftypes(entry.message_update_rule)
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
        entry.message_update_rule = first(applicable_rules)
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
    entry_recognition_factor = recognitionFactor(entry.interface.edge) # Recognition factor for outbound edge
    recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in entry.interface.node.interfaces
        node_interface_recognition_factor = recognitionFactor(node_interface.edge)

        if node_interface === entry.interface
            push!(inbound_types, Nothing)
        elseif node_interface_recognition_factor === entry_recognition_factor
            # Edge is internal, accept message
            push!(inbound_types, inferred_outbound_types[node_interface.partner])
        elseif !(node_interface_recognition_factor in recognition_factors)
            # Edge is external, accept marginal (if marginal is not already accepted)
            push!(inbound_types, ProbabilityDistribution) 
        end

        push!(recognition_factors, node_interface_recognition_factor)
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
    input_type_validators = 
        String["length(input_types) == $(length(inbound_types.args))"]
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

"""
Construct argument code for structured VB updates
"""
collectInbounds(entry::ScheduleEntry, ::Type{T}) where T<:StructuredVariationalRule = collectStructuredVariationalNodeInbounds(entry.interface.node, entry)

"""
Construct the inbound code that computes the message for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectStructuredVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry)
    interface_to_schedule_entry = current_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_algorithm.target_to_marginal_entry

    inbounds = Any[]
    entry_recognition_factor = recognitionFactor(entry.interface.edge)
    local_clusters = localRecognitionFactorization(entry.interface.node)

    recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        node_interface_recognition_factor = recognitionFactor(node_interface.edge)

        if node_interface === entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, nothing)
        elseif (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, ProbabilityDistribution))
        elseif node_interface_recognition_factor === entry_recognition_factor
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif !(node_interface_recognition_factor in recognition_factors)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            target = local_clusters[node_interface_recognition_factor]
            push!(inbounds, target_to_marginal_entry[target])
        end

        push!(recognition_factors, node_interface_recognition_factor)
    end

    return inbounds
end