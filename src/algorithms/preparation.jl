###########################################
# Shared methods for algorithm construction
###########################################

parameters{T<:ProbabilityDistribution}(message_type::Type{Message{T}}) = message_type.parameters[1].parameters

parameters{T<:ProbabilityDistribution}(distribution_type::Type{T}) = distribution_type.parameters

parameters(type_var::TypeVar) = type_var.ub.parameters

function parameters(data_type::DataType)
    # Extract parameters of type ForneyLab.Message{T<:ForneyLab.MvGaussianDistribution{dims}}
    # TODO: find more specific calling signature
    if data_type <: Message
        return data_type.parameters[1].ub.parameters
    else
        error("parameters function not valid for $(data_type)")
    end
end

function extractParameters(method_signature::SimpleVector, call_signature::Vector{DataType})
    # Constructs a dictionary of parameter variables in method_signature mapped to values in call_signature

    params_dict = Dict()
    for p in 4:length(method_signature) # Loop of indices of inbound types
        if call_signature[p] != Void # Skip irrelevant inbounds
            method_params = parameters(method_signature[p])
            call_params = parameters(call_signature[p])
            for r in 1:length(method_params) # Iterate over all found parameters for that specific inbound and the pairs to the dictionary
                push!(params_dict, method_params[r] => call_params[r])
            end
        end
    end

    return params_dict
end

function paramify(arr::Vector{Any})
    # Converts arr to a nice string for parameterized type construction
    str = ""
    for param in arr
        str *= "$(param), "
    end
    return str[1:end-2] # Chop off last comma
end

function extractOutboundType(outbound_arg)
    # Acceps the update rule argument of the outbound and returns a DataType with the parameterized outbound distribution type
    if typeof(outbound_arg) == TypeVar
        return outbound_arg.ub
    elseif typeof(outbound_arg) == DataType
        return outbound_arg
    else
        error("outbound_arg of unrecognized type: $(typeof(outbound_arg))")
    end
end

function collectAllOutboundTypes(rule::Function, call_signature::Vector, node::Node)
    # Find the available methods for the update function 'rule' that satisfy 'call_signature' and push the result to 'outbound_types'
    # Note the node::Node argument, so for specific node types this function may be overloaded

    outbound_types = DataType[]

    available_methods = methods(rule, call_signature)
    for method in available_methods

        method_signature = method.sig.types
        outbound_type = extractOutboundType(method_signature[3]) # Third entry is always the outbound distribution

        if typeof(node) <: TerminalNode
            push!(outbound_types, typeof(node.value))
        elseif isempty(parameters(outbound_type)) # Are there no type variables defined for the outbound type?
            push!(outbound_types, outbound_type) # Simply push the found outbound type on the stack
        else  # Outbound type has parameters that need to be inferred
            params_dict = extractParameters(method_signature, call_signature) # Extract parameters and values of inbound types
            outbound_params = parameters(outbound_type) # Extract parameters of outbound type

            # Construct new outbound type definition with substituted values
            param_values = Any[]
            for param in outbound_params
                push!(param_values, params_dict[param])
            end

            substituted_outbound_type = eval(parse("$(outbound_type.name){$(paramify(param_values))}")) # Construct the type definition of the substituted outbund type
            push!(outbound_types, substituted_outbound_type) 
        end
    end

    return outbound_types
end

function inferOutboundType!(entry::ScheduleEntry)
    # Infers the outbound type from node and all available information on inbounds and post-processing
    inbound_types = entry.inbound_types
    outbound_interface_id = entry.outbound_interface_id
    node = entry.node

    # Find all outbound types compatible with entry.rule
    outbound_types = collectAllOutboundTypes(entry.rule, [typeof(node); Type{Val{outbound_interface_id}}; Any; inbound_types], node)

    # The outbound outbound_types should contain just one element (there should be just one available update rule)
    if length(outbound_types) == 0
        error("No calculation rule available for inbound types $(inbound_types) on node $(node)")
    elseif length(outbound_types) > 1
        error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types) on node $(node)")
    elseif outbound_types[1] <: ProbabilityDistribution
        # There is only one possible outbound type and it is a probability distribution
        entry.intermediate_outbound_type = outbound_types[1]
    else
        error("Unknown output of message calculation rule: $(outbound_types[1]) for node $(node)")
    end

    
    if isdefined(entry, :post_processing)
        entry.outbound_type = inferOutboundTypeAfterPostProcessing(entry) # If post-processing is defined, entry.outbound_type might differ from entry.intermediate_outbound_type
    else
        entry.outbound_type = entry.intermediate_outbound_type # No post-processing, outbound type remains unaltered
    end

    return entry
end

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

function buildExecute!(entry::ScheduleEntry, inbound_arguments::Vector)
    # Constructs the execute function with optional post processing folded in
    # Additionally, this function completes the rule arguments for the update function call with the pointer to the outbound distribution

    if !isdefined(entry, :post_processing)
        # Create pointer to outbound distribution
        outbound_dist = entry.node.interfaces[entry.outbound_interface_id].message.payload
        # No post-processing; assign the "compiled" computation rule as an anomynous function to entry.execute
        entry.execute = ( () -> entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments...) )
    else
        # Fold the post-processing operation into entry.execute()
        # Note that the distribution type after execute() in general does not correspond with the distribution type after post_processing().
        # Therefore, we need to copy the original disribution and perform the execute() and post_processing() operations on the duplicate.
        # Then we repopulate the fields of the original distribution with the fields of the duplicate.
        # This approach ensures that the distribution pointers remain valid.

        outbound_dist = entry.intermediate_outbound_type() # Create a dummy distribution for sumProductRule! to fill (is accessed from the entry.execute closure)

        entry.execute = ( () -> ( # Anonymous function for message passing and post-processing
            intermediate_outbound_distribution = entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments...); # In-place operation on previously created dummy distribution
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
