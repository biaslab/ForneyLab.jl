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
        error("No applicable msg update rule for $(cluster) with inbound types $(inbound_types)")
    elseif length(applicable_rules) > 1
        error("Multiple applicable msg update rules for $(cluster) with inbound types $(inbound_types): $(applicable_rules)")
    else
        marginal_update_rule = first(applicable_rules)
    end

    return marginal_update_rule
end

"""
Construct the marginal computations table for a given recognition factor
"""
function marginalTable(rf::RecognitionFactor)
    # Construct outbound types dictionary
    outbound_types = Dict{Interface, Type}()
    for entry in rf.schedule
        outbound_types[entry.interface] = outboundType(entry.message_update_rule)
    end

    variable_table = marginalTable(sort(collect(rf.variables)))
    cluster_table = [MarginalEntry(cluster, outbound_types) for cluster in sort(collect(rf.clusters))]

    return [variable_table; cluster_table]
end

"""
Find the inbound types that are required to compute a joint marginal over `target`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(cluster::Cluster, outbound_types::Dict{Interface, Type})
    inbound_types = Type[]
    cluster_recognition_factor = recognitionFactor(first(cluster.edges)) # Recognition factor for cluster
    recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in cluster.node.interfaces
        node_interface_recognition_factor = recognitionFactor(node_interface.edge) # Note: edges that are not assigned to a recognition factorization are assumed mean-field 

        if node_interface_recognition_factor === cluster_recognition_factor
            # Edge is internal, accept message
            push!(inbound_types, outbound_types[node_interface.partner])
        elseif !(node_interface_recognition_factor in recognition_factors)
            # Edge is external, accept marginal (if marginal is not already accepted)
            push!(inbound_types, ProbabilityDistribution) 
        end

        push!(recognition_factors, node_interface_recognition_factor)
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
    input_type_validators = 
        String["length(input_types) == $(length(inbound_types.args))"]
    for (i, i_type) in enumerate(inbound_types.args)
        if i_type != :Nothing
            # Only validate inbounds required for update
            push!(input_type_validators, "ForneyLab.matches(input_types[$i], $i_type)")
        end
    end

    expr = parse("""
        begin
            mutable struct $name <: MarginalRule{$node_type} end
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{<:Type}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end