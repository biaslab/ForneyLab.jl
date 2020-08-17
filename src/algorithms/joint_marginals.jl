export
MarginalRule,
@marginalRule

"""
`MarginalRule{factor_type}` specifies a joint marginal update rule with respect
to a node of type `factor_type`.
"""
abstract type MarginalRule{factor_type} <: MarginalUpdateRule end

"""
Construct a `MarginalEntry` for computing the marginal over `cluster`
through a node-specific joint marginal update rule.
"""
function MarginalEntry(target::Cluster, outbound_types::Dict{Interface, Type})
    inbound_types = collectInboundTypes(target, outbound_types)
    marginal_update_rule = inferMarginalRule(target, inbound_types)

    # Collect inbound interfaces
    inbound_interfaces = Interface[]
    for edge in target.edges
        if edge.a in target.node.interfaces
            push!(inbound_interfaces, edge.a.partner) # Partner is the required inbound interface
        else
            push!(inbound_interfaces, edge.b.partner)
        end
    end

    return MarginalEntry(target, inbound_interfaces, marginal_update_rule)
end

"""
Infer the rule that computes the joint marginal over `cluster`
"""
function inferMarginalRule(cluster::Cluster, inbound_types::Vector{<:Type})
    # Find applicable rule(s)
    applicable_rules = Type[]
    for rule in leaftypes(MarginalRule{typeof(cluster.node)})
        if isApplicable(rule, inbound_types)
            push!(applicable_rules, rule)
        end
    end

    # Select and set applicable rule
    if isempty(applicable_rules)
        error("No applicable marginal update rule for $(typeof(cluster.node)) node with inbound types: $(join(inbound_types, ", "))")
    elseif length(applicable_rules) > 1
        error("Multiple applicable marginal update rules for $(typeof(cluster.node)) node with inbound types: $(join(inbound_types, ", "))")
    else
        marginal_update_rule = first(applicable_rules)
    end

    return marginal_update_rule
end

"""
Construct the marginal computations table for a given posterior factor.
The marginal table defines which marginal posteriors are computed.
"""
function marginalTable(pf::PosteriorFactor)
    # Construct outbound types dictionary
    outbound_types = Dict{Interface, Type}()
    for entry in pf.schedule
        outbound_types[entry.interface] = outboundType(entry.message_update_rule)
    end

    variable_table = marginalTable(sort(collect(pf.target_variables)))
    cluster_table = [MarginalEntry(cluster, outbound_types) for cluster in sort(collect(pf.target_clusters))]

    return [variable_table; cluster_table]
end

"""
Find the inbound types that are required to compute a joint marginal over `target`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(cluster::Cluster, outbound_types::Dict{Interface, Type})
    inbound_types = Type[]
    cluster_pf = posteriorFactor(first(cluster.edges))
    encountered_external_regions = Set{Region}()
    for node_interface in cluster.node.interfaces
        current_region = region(cluster.node, node_interface.edge) # Returns a Variable if no cluster is assigned
        current_pf = posteriorFactor(node_interface.edge) # Returns an Edge if no posterior factor is assigned
        inbound_interface = ultimatePartner(node_interface)

        if (current_pf === cluster_pf) && (current_region === cluster)
            # Edge is internal and in cluster, accept message
            push!(inbound_types, outbound_types[inbound_interface])
        elseif (current_pf === cluster_pf)
            # Edge is internal but not in cluster, signal to rule signature that edge will be marginalized out by appending "Nothing"
            push!(inbound_types, Nothing)
        elseif !(current_region in encountered_external_regions)
            # Edge is external and cluster is not yet encountered, accept probability distribution
            push!(inbound_types, ProbabilityDistribution)
            push!(encountered_external_regions, current_region) # Register current region with encountered external regions
        end
    end

    return inbound_types
end

"""
`@marginalRule` registers a marginal update rule for a (joint) marginal
by defining the rule type and the corresponding methods for the `isApplicable` function.
If no name (type) for the new rule is passed, a unique name (type) will be generated.
Returns the rule type.
"""
macro marginalRule(fields...)
    # Init required fields in macro scope
    node_type = :unknown
    inbound_types = :unknown
    name = :auto # Triggers automatic naming unless overwritten

    # Loop over fields because order is unknown
    for arg in fields
        (arg.args[1] == :(=>)) || error("Invalid call to @marginalRule")

        if arg.args[2].value == :node_type
            node_type = arg.args[3]
        elseif arg.args[2].value == :inbound_types
            inbound_types = arg.args[3]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[2].value == :name
            name = arg.args[3]
        else
            error("Unrecognized field $(arg.args[2].value) in call to @marginalRule")
        end
    end

    # Assign unique name if not set already
    if name == :auto
        # Added hash ensures that the rule name is unique
        msg_types_hash = string(hash(vcat([outbound_type], inbound_types)))[1:6]
        name = Symbol("Marginal$(node_type)$(msg_types_hash)")
    end

    # Build validators for isApplicable
    input_type_validators = Expr[]

    push!(input_type_validators, :(length(input_types) == $(length(inbound_types.args))))
    for (i, i_type) in enumerate(inbound_types.args)
        if i_type != :Nothing
            # Only validate inbounds required for update
            push!(input_type_validators, :(ForneyLab.matches(input_types[$i], $i_type)))
        end
    end

    expr = quote
        struct $name <: MarginalRule{$node_type} end
        ForneyLab.isApplicable(::Type{$name}, input_types::Vector{<:Type}) = begin
            $(reduce((current, item) -> :($current && $item), input_type_validators, init = :true))
        end
    end

    return esc(expr)
end

"""
Construct the inbound code that computes the marginal for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
collectInbounds(entry::MarginalEntry) = collectMarginalNodeInbounds(entry.target.node, entry)

function collectMarginalNodeInbounds(::FactorNode, entry::MarginalEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry
    inbound_cluster = entry.target # Entry target is a cluster

    inbounds = Any[]
    entry_pf = posteriorFactor(first(entry.target.edges))
    encountered_external_regions = Region[]
    for node_interface in entry.target.node.interfaces
        current_region = region(inbound_cluster.node, node_interface.edge) # Note: edges that are not assigned to a posterior factor are assumed mean-field
        current_pf = posteriorFactor(node_interface.edge) # Returns an Edge if no posterior factor is assigned
        inbound_interface = ultimatePartner(node_interface)

        if isClamped(inbound_interface)
            # Edge is clamped, hard-code marginal of constant node
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif (current_pf === entry_pf)
            # Edge is internal, collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif !(current_region in encountered_external_regions)
            # Edge is external and region is not yet encountered, collect marginal from marginal dictionary
            push!(inbounds, target_to_marginal_entry[current_region])
            push!(encountered_external_regions, current_region) # Register current region with encountered external regions
        end
    end

    return inbounds
end
