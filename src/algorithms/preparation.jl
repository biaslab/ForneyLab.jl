###########################################
# Shared methods for algorithm construction
###########################################

function inferOutboundType!(entry::ScheduleEntry, node::Node, allowed_rules::Vector{Function})
    # Infers the outbound type from node and all available information on inbounds and post-processing
    inbound_types = entry.inbound_types
    outbound_interface_id = entry.outbound_interface_id

    # Find all compatible calculation rules for the SumProduct algorithm
    outbound_types = []
    for update_function in allowed_rules
        available_rules = methods(update_function, [typeof(node); Type{Val{outbound_interface_id}}; inbound_types; Any])
        for rule in available_rules
            push!(outbound_types, rule.sig.types[end]) # Push the found outbound type on the stack
        end
    end

    # The outbound outbound_types should contain just one element (there should be just one available update rule)
    if length(outbound_types) == 0
        error("No calculation rule available for inbound types $(inbound_types) on node $(node)")
    elseif length(outbound_types) > 1
        error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types) on node $(node)")
    elseif outbound_types[1] == Any
        # The computation rule produces Any, which indicates that the node is parametrized by its outbound type (e.g. a TerminalNode)
        (typeof(node).parameters[1] <: ProbabilityDistribution) || error("$(typeof(node)) $(node.id) must be parametrized by a ProbabilityDistribution")
        entry.intermediate_outbound_type = typeof(node).parameters[1]
    elseif outbound_types[1] <: ProbabilityDistribution
        # There is only one possible outbound type and it is a probability distribution
        entry.intermediate_outbound_type = outbound_types[1]
    else
        error("Unknown output of message calculation rule: $(outbound_types[1]) for node $(node)")
    end

    
    if isdefined(entry, :post_processing)
        entry.outbound_type = inferOutboundTypeAfterPostProcessing(entry) # If post-processing is defined, entry.outbound_type might differ from entry.intermediate_outbound_type
    else
        entry.outbound_type = entry.intermediate_outbound_type # No post-processing
    end

    return entry
end

inferOutboundType!(entry::ScheduleEntry, allowed_rules::Vector{Function}) = inferOutboundType!(entry, entry.node, allowed_rules)

function inferOutboundTypeAfterPostProcessing(entry::ScheduleEntry)
    outbound_types = []
    available_post_processing_rules = methods(entry.post_processing, [Any, entry.intermediate_outbound_type])

    length(available_post_processing_rules) > 0 || error("No post processing available as $(entry.post_processing) on $(entry.intermediate_outbound_type). Please make sure all post-processing functions fit the signature: function{T}(::Type{T<:ProbabilityDistribution}, d::ProbabilityDistribution)")

    for post_processing_rule in available_post_processing_rules
        outbound_type = post_processing_rule.sig.types[1].parameters[1]
        (outbound_type <: ProbabilityDistribution) || continue # Skip when result is not a distribution
        push!(outbound_types, outbound_type) # Push the found outbound type (first entry) on the stack
    end

    length(outbound_types) <= 1 || error("Multiple post processing possibilities for $(entry.post_processing): ($(outbound_types)).")

    return outbound_types[1]
end


#######################################################
# Shared methods for algorithm preparation/compilation
#######################################################

function buildExecute!(entry::ScheduleEntry, rule_arguments::Vector)
    # Constructs the execute function with optional post processing folded in
    # Additionally, this function completes the rule arguments for the update function call with the pointer to the outbound distribution

    if !isdefined(entry, :post_processing)
        # Add outbound distribution to rule_arguments
        push!(rule_arguments, entry.node.interfaces[entry.outbound_interface_id].message.payload)
        # No post-processing; assign the "compiled" computation rule as an anomynous function to entry.execute
        entry.execute = ( () -> entry.rule(entry.node, Val{entry.outbound_interface_id}, rule_arguments...) )
    else
        # Fold the post-processing operation into entry.execute()
        # Note that the distribution type after execute() in general does not correspond with the distribution type after post_processing().
        # Therefore, we need to copy the original disribution and perform the execute() and post_processing() operations on the duplicate.
        # Then we repopulate the fields of the original distribution with the fields of the duplicate.
        # This approach ensures that the distribution pointers remain valid.

        push!(rule_arguments, entry.intermediate_outbound_type()); # Append a dummy distribution for sumProduct! to fill

        entry.execute = ( () -> ( # Anonymous function for message passing and post-processing
            intermediate_outbound_distribution = entry.rule(entry.node, Val{entry.outbound_interface_id}, rule_arguments...); # In-place operation on previously created dummy distribution
            new_outbound_distribution = entry.post_processing(entry.outbound_type, intermediate_outbound_distribution); # Not an in-place operation
            original_outbound_distribution = entry.node.interfaces[entry.outbound_interface_id].message.payload; # Get original pointer to outbound distribution on interface
            
            # Duplicate parameters of new into original
            injectParameters!(original_outbound_distribution, new_outbound_distribution);
            
            return original_outbound_distribution
        ))
    end

    return entry
end

function injectParameters!{T<:ProbabilityDistribution}(destination::T, source::T)
    # Fill the parameters of a destination distribution with the copied parameters of the source

    for field in fieldnames(source)
        setfield!(destination, field, deepcopy(getfield(source, field)))
    end

    return destination
end