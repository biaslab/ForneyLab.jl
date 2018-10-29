export
MarginalRule,
@marginalRule

"""
MarginalRule{factor_type} specifies a joint marginal update rule with respect
to a node of type `factor_type`.
"""
abstract type MarginalRule{factor_type} <: MarginalUpdateRule end

"""
Construct a MarginalScheduleEntry for computing the marginal over `cluster`
through a node-specific joint marginal update rule.
"""
function MarginalScheduleEntry(cluster::Cluster, outbound_types::Dict{Interface, Type})
    inbound_types = collectInboundTypes(cluster, outbound_types)
    marginal_update_rule = inferMarginalRule(cluster, inbound_types)
    
    # Collect inbound interfaces 
    inbound_interfaces = Interface[]
    for edge in cluster.edges
        if edge.a in cluster.node.interfaces
            push!(inbound_interfaces, edge.a.partner) # Partner is the required inbound interface
        else
            push!(inbound_interfaces, edge.b.partner)
        end
    end

    return MarginalScheduleEntry(cluster, inbound_interfaces, marginal_update_rule)
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
Find the inbound types that are required to compute a joint marginal over `target`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(cluster::Cluster, outbound_types::Dict{Interface, Type})
    inbound_types = Type[]
    cluster_recognition_factor_id = recognitionFactorId(first(cluster.edges)) # Recognition factor id for cluster
    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in cluster.node.interfaces
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge) # Note: edges that are not assigned to a recognition factorization are assumed mean-field 

        if node_interface_recognition_factor_id == cluster_recognition_factor_id
            # Edge is internal, accept message
            push!(inbound_types, outbound_types[node_interface.partner])
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Edge is external, accept marginal (if marginal is not already accepted)
            push!(inbound_types, ProbabilityDistribution) 
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbound_types
end

function marginalSchedule(q_factors::Vector{RecognitionFactor}, schedule::Schedule)
    # Construct outbound types dictionary
    outbound_types = Dict{Interface, Type}()
    for entry in schedule
        outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    # Construct marginal schedule
    marginal_schedule = MarginalScheduleEntry[]
    for q_factor in q_factors
        # Construct schedule for computing marginals over variables
        variable_schedule = [MarginalScheduleEntry(variable) for variable in sort(collect(q_factor.variables))]
        marginal_schedule = [marginal_schedule; variable_schedule]

        # Construct schedule for computing marginals over clusters
        cluster_schedule = [MarginalScheduleEntry(cluster, outbound_types) for cluster in sort(collect(q_factor.clusters))]
        marginal_schedule = [marginal_schedule; cluster_schedule]
    end

    return marginal_schedule
end

marginalSchedule(q_factor::RecognitionFactor, schedule::Schedule) = marginalSchedule([q_factor], schedule)

"""
@marginalRule registers a marginal update rule for a (joint) marginal
by defining the rule type and the corresponding methods for the isApplicable functions.
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
    input_type_validators = String[]
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