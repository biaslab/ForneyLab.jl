###########################################
# Shared methods for algorithm construction
###########################################

parameters{T<:ProbabilityDistribution}(message_type::Type{Message{T}}) = message_type.parameters[1].parameters

parameters{T<:ProbabilityDistribution}(distribution_type::Type{T}) = distribution_type.parameters

parameters(type_var::TypeVar) = type_var.ub.parameters

function parameters(data_type::DataType)
    # Extract parameters of type ForneyLab.Message{T<:ForneyLab.MvGaussianDistribution{dims}}
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

paramify(arr::Vector{Any}) = join(arr, ", ")

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

            parametrized_outbound_type = eval(parse("$(outbound_type.name){" * join(param_values,", ") * "}")) # Construct the parametrized outbund type
            push!(outbound_types, parametrized_outbound_type)
        end
    end

    return outbound_types
end

function inferOutboundType!(entry::ScheduleEntry)
    # Infer the outbound type from the node type and the types of the inbounds

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
        entry.outbound_type = outbound_types[1]
        entry.rule_is_approximate = false
    else
        error("Unknown output of message calculation rule: $(outbound_types[1]) for node $(node)")
    end

    return entry
end


#######################################################
# Shared methods for algorithm preparation/compilation
#######################################################

function buildExecute!(entry::ScheduleEntry, inbound_arguments::Vector)
    # Construct the entry.execute function.
    # This function is called by the prepare methods of inference algorithms.

    # Get pointer to the outbound distribution
    outbound_dist = entry.node.interfaces[entry.outbound_interface_id].message.payload

    # Save the "compiled" message computation rule as an anomynous function in entry.execute
    entry.execute = ( () -> entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments...) )

    return entry
end

function injectParameters!{T<:ProbabilityDistribution}(destination::T, source::T)
    # Fill the parameters of a destination distribution with the copied parameters of the source

    for field in fieldnames(source)
        setfield!(destination, field, deepcopy(getfield(source, field)))
    end

    return destination
end
