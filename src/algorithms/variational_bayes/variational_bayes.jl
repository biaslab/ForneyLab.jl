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
    # TODO: incorporate structured updates

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
                                                rule_type::Type{T},
                                                ::Dict{Interface, DataType})
    # Find outbound id
    outbound_id = findfirst(entry.interface.node.interfaces, entry.interface)    
    
    # Find applicable rule(s)
    applicable_rules = DataType[]
    for rule in leaftypes(entry.msg_update_rule)
        if isApplicable(rule, outbound_id)
            push!(applicable_rules, rule)
        end
    end

    # Select and set applicable rule
    if isempty(applicable_rules)
        error("No applicable msg update rule for $(entry) with outbound id $(outbound_id)")
    elseif length(applicable_rules) > 1
        error("Multiple applicable msg update rules for $(entry) with outbound id $(outbound_id)")
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
    outbound_id = :unknown # Mean-field variational rule does not depend on inbounds; only on node type and outbound id
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.head == :(=>)) || error("Invalid call to @variationalRule")

        if arg.args[1].args[1] == :node_type
            node_type = arg.args[2]
        elseif arg.args[1].args[1] == :outbound_type
            outbound_type = arg.args[2]
            (outbound_type.head == :curly && outbound_type.args[1] == :Message) || error("Outbound type for VariationalRule should be a Message")
        elseif arg.args[1].args[1] == :outbound_id
            outbound_id = arg.args[2]
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

    expr = parse("""
        begin
            type $name <: VariationalRule{$node_type} end
            ForneyLab.outboundType(::Type{$name}) = $outbound_type
            ForneyLab.isApplicable(::Type{$name}, outbound_id::Int64) = (outbound_id == $outbound_id)
            $name
        end
    """)

    return esc(expr)
end