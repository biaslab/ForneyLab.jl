export
NaiveVariationalRule,
variationalSchedule,
@naiveVariationalRule

abstract NaiveVariationalRule{factor_type} <: MessageUpdateRule

"""
variationalSchedule() generates a variational message passing schedule that computes the
marginals for each of the recognition distributions in the recognition factor.
"""
function variationalSchedule(recognition_factors::Vector{RecognitionFactor})
    # Schedule messages towards recognition distributions, limited to the internal edges
    schedule = ScheduleEntry[]
    nodes_connected_to_external_edges = Set{FactorNode}()
    for recognition_factor in recognition_factors
        schedule = [schedule; summaryPropagationSchedule(sort(collect(recognition_factor.variables)), limit_set=recognition_factor.internal_edges)]
        union!(nodes_connected_to_external_edges, nodesConnectedToExternalEdges(recognition_factor))
    end

    for entry in schedule
        if entry.interface.node in nodes_connected_to_external_edges
            local_recognition_factor_ids = localRecognitionFactorIds(entry.interface.node)
            if allunique(local_recognition_factor_ids) # Local recognition factorization is naive
                entry.msg_update_rule = NaiveVariationalRule{typeof(entry.interface.node)}
            else
                entry.msg_update_rule = StructuredVariationalRule{typeof(entry.interface.node)}
            end        
        else
            entry.msg_update_rule = SumProductRule{typeof(entry.interface.node)}
        end
    end

    inferUpdateRules!(schedule)

    return schedule
end

variationalSchedule(recognition_factor::RecognitionFactor) = variationalSchedule([recognition_factor])

function inferUpdateRule!{T<:NaiveVariationalRule}( entry::ScheduleEntry,
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
        error("Multiple applicable msg update rules for $(entry) with inbound types $(inbound_types)")
    else
        entry.msg_update_rule = first(applicable_rules)
    end

    return entry
end

function collectInboundTypes{T<:NaiveVariationalRule}(  entry::ScheduleEntry,
                                                        ::Type{T},
                                                        inferred_outbound_types::Dict{Interface, DataType})
    inbound_types = DataType[]
    for node_interface in entry.interface.node.interfaces
        if node_interface == entry.interface
            push!(inbound_types, Void)
        else
            # Edge is external, accept marginal
            push!(inbound_types, ProbabilityDistribution) 
        end
    end

    return inbound_types
end

"""
@naiveVariationalRule registers a variational update rule for the naive (mean-field)
factorization by defining the rule type and the corresponding methods for the 
outboundType and isApplicable functions. If no name (type) for the new rule is
passed, a unique name (type) will be generated. Returns the rule type.
"""
macro naiveVariationalRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.head == :(=>)) || error("Invalid call to @naiveVariationalRule")

        if arg.args[1].args[1] == :node_type
            node_type = arg.args[2]
        elseif arg.args[1].args[1] == :outbound_type
            outbound_type = arg.args[2]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for NaiveVariationalRule should be a Message")
        elseif arg.args[1].args[1] == :inbound_types
            inbound_types = arg.args[2]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[1].args[1] == :name
            name = arg.args[2]
        else
            error("Unrecognized field $(arg.args[1].args[1]) in call to @naiveVariationalRule")
        end
    end

    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("VB$(node_type)$(msg_types_hash)")
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
            type $name <: NaiveVariationalRule{$node_type} end
            ForneyLab.outboundType(::Type{$name}) = $outbound_type
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{DataType}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end