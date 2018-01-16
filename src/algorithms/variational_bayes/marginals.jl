export
MarginalUpdateRule,
MarginalScheduleEntry,
MarginalSchedule
@marginalRule

"""
A Cluster specifies a collection of `edges` adjacent to `node` that belong to the same
RecognitionFactor. Over a Cluster, a joint marginal can be computed.
"""
type Cluster
    id::Symbol
    node::FactorNode
    edges::Vector{Edge}

    function Cluster(node::FactorNode, edges::Vector{Edge})
        id = Symbol(join([edge.variable.id for edge in edges], "_"))
        self = new(id, node, edges)
        return self
    end
end

"""
A MarginalUpdateRule specifies how a marginal is calculated from incoming messages,
and an optional node function.
"""
abstract MarginalUpdateRule
abstract MarginalRule{factor_type} <: MarginalUpdateRule
abstract Product <: MarginalUpdateRule

"""
A `MarginalScheduleEntry` defines a marginal computation.
The `marginal_update_rule <: MarginalUpdateRule` defines the rule that is used
to calculate the (joint) marginal over `target`.
"""
type MarginalScheduleEntry
    target::Union{Cluster, Variable}
    interfaces::Vector{Interface}
    marginal_update_rule::DataType
end

typealias MarginalSchedule Vector{MarginalScheduleEntry}

"""
marginalSchedule() generates a marginal schedule that computes the
marginals for each target argument.
"""
function marginalSchedule(targets::Vector{Union{Variable, Cluster}})
    marginal_schedule = MarginalScheduleEntry[]

    # Schedule a marginal update rule for each of the targets
    for target in targets
        if isa(target, Variable)
            # Target is over a single variable, compute marginal through multiplication
            target_edge = first(target.edges) # For the sake of consistency, we always take the first edge.
            if target_edge.a == nothing # First handle cases where there is a `dangling` edge
                entry = MarginalScheduleEntry(target, [target_edge.b], Void)
            elseif target_edge.b == nothing
                entry = MarginalScheduleEntry(target, [target_edge.a], Void)
            else
                entry = MarginalScheduleEntry(target, [target_edge.a, target_edge.b], Product)
            end
        else
            # Target is over multiple edges, compute joint marginal through specific marginal update rule
            inbound_types = collectInboundTypes(target, inferred_outbound_types) # TODO: how to obtain inferred_outbound_types
            marginal_update_rule = inferMarginalRule(target, inbound_types)
            
            # Collect inbound interfaces 
            target_interfaces = Interface[]
            for target_edge in target.edges
                if target_edge.a in target.node.interfaces
                    push!(target_interfaces, target_edge.b) # Push partner to target interfaces
                else
                    push!(target_interfaces, target_edge.a)
                end
            end

            entry = MarginalScheduleEntry(target, target_interfaces, marginal_update_rule)
        end

        push!(marginal_schedule, entry)
    end

    return marginal_schedule
end

function inferMarginalRule(cluster::Cluster, inbound_types::Vector{DataType})
    # Find applicable rule(s)
    applicable_rules = DataType[]
    for rule in leaftypes(MarginalRule{typeof(cluster.node)})
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
        marginal_update_rule = first(applicable_rules)
    end

    return marginal_update_rule
end

"""
Find the inbound types that are required to compute a joint marginal over `target`.
Returns a vector with inbound types that correspond with required interfaces.
"""
function collectInboundTypes(cluster::Cluster, inferred_outbound_types::Dict{Interface, DataType})
    inbound_types = DataType[]
    cluster_recognition_factor_id = recognitionFactorId(first(cluster.edges)) # Recognition factor id for cluster
    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in cluster.node.interfaces
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge) # Note: edges that are not assigned to a recognition factorization are assumed mean-field 

        if node_interface_recognition_factor_id == cluster_recognition_factor_id
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
        (arg.head == :(=>)) || error("Invalid call to @marginalRule")

        if arg.args[1].args[1] == :node_type
            node_type = arg.args[2]
        elseif arg.args[1].args[1] == :inbound_types
            inbound_types = arg.args[2]
            (inbound_types.head == :tuple) || error("Inbound types should be passed as Tuple")
        elseif arg.args[1].args[1] == :name
            name = arg.args[2]
        else
            error("Unrecognized field $(arg.args[1].args[1]) in call to @marginalRule")
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
        if i_type != :Void
            # Only validate inbounds required for update
            push!(input_type_validators, "ForneyLab.matches(input_types[$i], $i_type)")
        end
    end

    expr = parse("""
        begin
            type $name <: MarginalRule{$node_type} end
            ForneyLab.isApplicable(::Type{$name}, input_types::Vector{DataType}) = $(join(input_type_validators, " && "))
            $name
        end
    """)

    return esc(expr)
end