export
VariationalRule,
variationalSchedule,
@variationalRule

abstract VariationalRule{factor_type} <: MessageUpdateRule

"""
variationalSchedule() generates a variational message passing schedule that computes the
marginals for each of the recognition distributions in the recognition factor.
"""
function variationalSchedule(recognition_factor::RecognitionFactor)
    # TODO: more efficient scheduling (see e.g. commit a5272c in master branch)

    internal_edges = recognition_factor.internal_edges
    # Schedule messages towards recognition distributions, limited to the internal edges
    schedule = summaryPropagationSchedule(sort(collect(recognition_factor.variables)), limit_set=internal_edges)

    # external_edges are the difference between all edges connected to nodes, and the internal edges
    subgraph_nodes = nodes(internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), internal_edges)
    # nodes_connected_to_external_edges are the nodes connected to external edges that are also connected to internal edges
    nodes_connected_to_external_edges = intersect(nodes(external_edges), subgraph_nodes)
    for entry in schedule
        if entry.interface.node in nodes_connected_to_external_edges
            entry.msg_update_rule = VariationalRule{typeof(entry.interface.node)}
        else
            entry.msg_update_rule = SumProductRule{typeof(entry.interface.node)}
        end
    end

    inferUpdateRules!(schedule)

    return schedule
end

function inferUpdateRule!{T<:VariationalRule}(  entry::ScheduleEntry,
                                                ::Type{T},
                                                inferred_outbound_types::Dict{Interface, DataType})
    # Collect inbound types
    inbound_types = DataType[]
    for node_interface in entry.interface.node.interfaces
        if is(node_interface, entry.interface)
            push!(inbound_types, Void)
        elseif is(recognitionFactor(node_interface.edge), recognitionFactor(entry.interface.edge))
            # Both edges are internal in the same recognition factor, require message
            push!(inbound_types, inferred_outbound_types[node_interface.partner])
        else
            # A recognition factor border is crossed, require recognition distribution
            # for the variable corresponding to the external edge
            push!(inbound_types, typeof(current_recognition_factorization.recognition_distributions[node_interface.edge.variable]))
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

"""
@variationalRule registers a variational update rule by defining the rule type
and the corresponding methods for the outboundType and isApplicable functions.
If no name (type) for the new rule is passed, a unique name (type) will be
generated. Returns the rule type.
"""
macro variationalRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    outbound_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.head == :(=>)) || error("Invalid call to @variationalRule")

        if arg.args[1].args[1] == :node_type
            node_type = arg.args[2]
        elseif arg.args[1].args[1] == :outbound_type
            outbound_type = arg.args[2]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for VariationalRule should be a Message")
        elseif arg.args[1].args[1] == :inbound_types
            inbound_types = arg.args[2]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[1].args[1] == :name
            name = arg.args[2]
        else
            error("Unrecognized field $(arg.args[1].args[1]) in call to @variationalRule")
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
            push!(input_type_validators, "(input_types[$i]==$i_type)")
        end
    end

    expr = parse("""
        begin
            type $name <: VariationalRule{$node_type} end
            ForneyLab.outboundType(::Type{$name}) = $outbound_type
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{DataType}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end