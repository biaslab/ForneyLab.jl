###########################################
# Shared methods for algorithm construction
###########################################

function collectTypeVarValues(method::Method, argument_types::Vector{DataType}, results::Dict{Symbol,Any}=Dict{Symbol,Any}())
    # Resolve all the values of TypeVar objects in method.sig.types, based on the types of the actual arguments (argument_types)
    # The results are collected in a Dict, indexed by the names of the TypeVars.
    # Example: if method.sig.types[2] = Message{MvGaussian{dims}}
    #           and argument_types[2] = Message{MvGaussian{3}}
    #          {:dims => 3} is added to the results Dict.

    for a=1:length(method.sig.types)
        _collectTypeVarValues!(method.sig.types[a], argument_types[a], results)
    end

    return results
end

function _collectTypeVarValues!(method_sig_type::TypeVar, param_value, results::Dict{Symbol,Any})
    # Internal function called by collectTypeVarValues
    results[method_sig_type.name] = param_value
end

function _collectTypeVarValues!(method_sig_type::DataType, arg_type::DataType, results::Dict{Symbol,Any})
    # Internal function called by collectTypeVarValues
    # The recursion performs a depth-first search through the (hierarchical) parameters of the argument type
    if length(method_sig_type.parameters) > 0
        for p = 1:length(method_sig_type.parameters)
            _collectTypeVarValues!(method_sig_type.parameters[p], arg_type.parameters[p], results)
        end
    end
end

function containsTypeVars(dtype::DataType)
    # Check if dtype contains TypeVar parameters
    for param in dtype.parameters
        if isa(param, TypeVar)
            return true
        elseif isa(param, DataType)
            containsTypeVars(param) && return true
        end
    end

    return false
end

function resolveTypeVars(dtype::TypeVar, method::Method, argument_types::Vector{DataType})
    # Find the value of dtype based on the method and its instantiation argument types
    # Return: (resolved_dtype, completely_resolved)
    typevar_values = collectTypeVarValues(method, argument_types)
    try
        return (typevar_values[dtype.name], true)
    catch
        (dtype, false)
    end
end

function resolveTypeVars(dtype::DataType, method::Method, argument_types::Vector{DataType})
    # Build a DataType based on dtype, but with all TypeVar parameters replaced by their corresponding values.
    # The values are determined from the method specification and the calling signature of the method.
    # Return: (resolved_dtype, completely_resolved)

    containsTypeVars(dtype) || return (dtype, true) # Nothing to resolve
    typevar_values = collectTypeVarValues(method, argument_types) # Collect all TypeVar values in a Dict indexed by the TypeVar name

    return resolveTypeVars(dtype, typevar_values) # TODO: implement this method to actually replace TypeVars with their corresponding values
end

function resolveTypeVars!(dtype::DataType, lookup_table::Dict{Symbol,Any})
    completely_resolved = true
    error("TODO")
end

function collectAllOutboundTypes(rule::Function, call_signature::Vector, node::Node)
    # Collect all outbound distribution types that can be generated with the specified rule and calling_signature combination
    # Note the node argument: this function can be overloaded for specific nodes
    # Returns: outbound_types::Vector{DataType}
    # The entries of outbound_types can be <: ProbabilityDistribution or Approximation

    outbound_types = DataType[]

    for method in methods(rule, call_signature)
        ob_type = method.sig.types[3] # Third entry is always the outbound distribution

        if typeof(node) <: TerminalNode
            ob_type = typeof(node.value)
        else
            substituteParameterValues(ob_type, method, call_signature)
        end

        # Is the method an approximate msg calculation rule?
        if (typeof(method.sig.types[end])==DataType
            && length(method.sig.types[end].parameters)==1
            && method.sig.types[end].parameters[1]<:ApproximationType)

            ob_type = Approximation{ob_type,method.sig.types[end].parameters[1]}
        end

        push!(outbound_types, ob_type)
    end

    return outbound_types
end

function collectAllOutboundTypes(rule::Function, call_signature::Vector, node::Union{GainNode, GainAdditionNode, GainEqualityNode})
    # Outbound type collection overloading for nodes with an (optional) fixed gain.
    # This is necessary because the fixed gain determines the dimensionality of the outbound distribution.

    outbound_types = DataType[]

    for method in methods(rule, call_signature)
        ob_type = method.sig.types[3] # Third entry is always the outbound distribution

        if !isempty(parameters(ob_type)) # Outbound type has parameters that need to be inferred
            params_dict = extractParameterValues(method, call_signature) # Extract parameters and values of inbound types
            outbound_params = parameters(ob_type) # Extract parameters of outbound type

            # Construct new outbound type definition with substituted values
            param_values = Any[]
            for param in outbound_params
                if haskey(params_dict, param)
                    push!(param_values, params_dict[param])
                else # The outbound param is not available in the inbound parameter dictionary; we need to infer it from the fixed gain matrix
                    if param.name == :dims_n
                        push!(param_values, size(node.gain, 1))
                    elseif param.name == :dims_m
                        push!(param_values, size(node.gain, 2))
                    else
                        error("For the gain node with fixed gain, the dimensionalities in the calling signature need to be encoded as {dims_n, dims_m}")
                    end
                end
            end

            ob_type = eval(parse("$(ob_type.name){" * join(param_values, ", ") * "}")) # Construct the type definition of the substituted outbound type
        end

        # Is the method an approximate msg calculation rule?
        if (typeof(method.sig.types[end])==DataType
            && length(method.sig.types[end].parameters)==1
            && method.sig.types[end].parameters[1]<:ApproximationType)

            ob_type = Approximation{ob_type,method.sig.types[end].parameters[1]}
        end

        push!(outbound_types, ob_type)
    end

    return outbound_types
end

function buildRuleSignature(rule::Function, node::Node, outbound_interface_id::Int64, inbound_types::Vector{DataType}, approx=false)
    if approx == false
        return vcat([typeof(node); Type{Val{outbound_interface_id}}; Any], inbound_types)
    else
        return vcat([typeof(node); Type{Val{outbound_interface_id}}; Any], inbound_types, approx)
    end
end

function collectOutboundTypes(rule::Function, node::Node, outbound_interface_id::Int64, inbound_types::Vector{DataType}, approx=false; unroll_partitioned_inbounds::Bool=false)
    # Collect all possible outbound types for a given rule, node, and inbound types
    if unroll_partitioned_inbounds
        num_factors = partitionedInboundTypesSize(inbound_types)
        if num_factors > 0
            # Try to find extra outbound types by splitting partitioned inbounds
            modified_inbound_types = stripPartitionedInboundTypes(inbound_types)
            sig = buildRuleSignature(rule, node, outbound_interface_id, modified_inbound_types, approx)
            return DataType[PartitionedDistribution{ob_type, num_factors} for ob_type in collectAllOutboundTypes(rule, sig, node)]
        else
            return DataType[]
        end
    else
        sig = buildRuleSignature(rule, node, outbound_interface_id, inbound_types, approx)
        return collectAllOutboundTypes(rule, sig, node)
    end
end

function inferOutboundType!(entry::ScheduleEntry)
    # Infer the outbound type from the node type and the types of the inbounds

    inbound_types = entry.inbound_types
    if isdefined(entry, :outbound_type)
        # The outbound type is already fixed, so we just need to validate that there exists a suitable computation rule

        # Try to find a matching exact rule
        for outbound_type in collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types)
            (entry.outbound_type <: outbound_type) && return(entry)
        end

        # Try to find a matching exact rule by unrolling partitioned inbounds
        for outbound_type in collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, unroll_partitioned_inbounds=true)
            if entry.outbound_type <: outbound_type
                entry.unrolling_factor = outbound_type.parameters[2]
                return(entry)
            end
        end

        # Try to find a matching approximate rule
        approx = isdefined(entry, :approximation) ? Type{entry.approximation} : Any
        for outbound_type in collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, approx)
            if entry.outbound_type <: outbound_type.parameters[1]
                entry.approximation = outbound_type.parameters[2]
                return entry
            end
        end

        # Try to find a matching approximate rule by unrolling partitioned inbounds
        for outbound_type in collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, approx, unroll_partitioned_inbounds=true)
            if entry.outbound_type <: outbound_type.parameters[1]
                entry.approximation = outbound_type.parameters[2]
                entry.unrolling_factor = outbound_type.parameters[1].parameters[2]
                return entry
            end
        end

        error("No suitable calculation rule available for schedule entry:\n$(entry)Inbound types: $(inbound_types).")
    else
        # Infer the outbound type

        # Try to apply an exact rule
        outbound_types = collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type.")
        end

        # Try to apply an exact rule by unrolling partitioned inbounds
        outbound_types = collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, unroll_partitioned_inbounds=true)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1]
            entry.unrolling_factor = outbound_types[1].parameters[2]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type.")
        end

        # Try to apply an approximate rule
        approx = isdefined(entry, :approximation) ? Type{entry.approximation} : Any
        outbound_types = collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, approx)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1].parameters[1]
            entry.approximation = outbound_types[1].parameters[2]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type and if required also an approximation method.")
        end

        # Try to apply an approximate rule by unrolling partitioned inbounds
        outbound_types = collectOutboundTypes(entry.rule, entry.node, entry.outbound_interface_id, inbound_types, approx, unroll_partitioned_inbounds=true)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1].parameters[1]
            entry.approximation = outbound_types[1].parameters[2]
            entry.unrolling_factor = outbound_types[1].parameters[1].parameters[2]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type and if required also an approximation method.")
        end

        error("No calculation rule available for schedule entry:\n$(entry)Inbound types: $(inbound_types).")
    end
end

function partitionedInboundTypesSize(inbound_types::Vector{DataType})
    # Return 0 if inbound_types contains no partitioned inbounds.
    # If a partitioned inbound is detected, return the number of factors.

    for ib_type in inbound_types
        if ib_type <: PartitionedDistribution
            return ib_type.parameters[2]
        elseif ((ib_type <: Message) && (ib_type.parameters[1] <: PartitionedDistribution))
            return ib_type.parameters[1].parameters[2]
        end
    end

    return 0
end

function stripPartitionedInboundTypes(inbound_types::Vector{DataType})
    # Simplify a list of inbound types that contains PartitionedDistribution or Message{PartitionedDistribution} entries.
    # Replace every partitioned inbound type by the type of the factors.
    # Check if all partitioned inbounds contain the same number of factors. If not, throw an error.
    stripped_inbound_types = Vector{DataType}()
    num_factors = 0
    for ib_type in inbound_types
        if ib_type <: PartitionedDistribution
            part_type = ib_type
        elseif ((ib_type <: Message) && (ib_type.parameters[1] <: PartitionedDistribution))
            part_type = ib_type.parameters[1]
        else
            # Regular inbound type
            push!(stripped_inbound_types, ib_type)
            continue
        end

        if num_factors == 0
            num_factors = part_type.parameters[2]
        elseif num_factors != part_type.parameters[2]
            error("PartitionedDistribution inbounds for a message computation rule should have the same number of factors.")
        end

        factor_type = part_type.parameters[1]
        push!(stripped_inbound_types, (ib_type <: Message) ? Message{factor_type} : factor_type)
    end

    return stripped_inbound_types
end

function interfacesFacingWrapsOrBuffers(graph::FactorGraph=currentGraph();
                                        include_wraps=true,
                                        include_buffers=true)
    # Return a list of interfaces in graph that face a wrap or writebuffer.
    interfaces = Interface[]

    # Collect wrap facing interfaces
    if include_wraps
        for wrap in wraps(graph)
            push!(interfaces, wrap.tail.interfaces[1].partner)
            if isdefined(graph, :block_size)
                push!(interfaces, wrap.head.interfaces[1].partner)
            end
        end
    end

    # Collect write buffer facing interfaces
    if include_buffers
        for entry in keys(graph.write_buffers)
            if typeof(entry) == Interface
                push!(interfaces, entry)
            elseif typeof(entry) == Edge
                push!(interfaces, entry.head)
                push!(interfaces, entry.tail)
            end
        end
    end

    return interfaces
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
    if entry.unrolling_factor > 0
        # Unroll partitioned inbounds and call message calculation rule for each factor
        factor_inbounds = copy(inbound_arguments) # inbounds for a single factor
        inbounds_to_unroll = Int64[]
        for i=1:length(inbound_arguments)
            if typeof(inbound_arguments[i]) <: PartitionedDistribution
                factor_inbounds[i] = inbound_arguments[i].factors[1]
                push!(inbounds_to_unroll, i)
            elseif ((typeof(inbound_arguments[i]) <: Message) && (typeof(inbound_arguments[i].payload) <: PartitionedDistribution))
                factor_inbounds[i] = Message(inbound_arguments[i].payload.factors[1])
                push!(inbounds_to_unroll, i)
            end
        end

        entry.execute = () -> begin
            for factor_idx=1:entry.unrolling_factor
                for i in inbounds_to_unroll
                    if typeof(factor_inbounds[i]) <: Message
                        factor_inbounds[i].payload = inbound_arguments[i].payload.factors[factor_idx]
                    else
                        factor_inbounds[i] = inbound_arguments[i].factors[factor_idx]
                    end
                end

                if isdefined(entry, :approximation)
                    entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist.factors[factor_idx], factor_inbounds..., entry.approximation)
                else
                    entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist.factors[factor_idx], factor_inbounds...)
                end
            end

            return outbound_dist
        end
    else
        if isdefined(entry, :approximation)
            entry.execute = ( () -> entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments..., entry.approximation) )
        else
            entry.execute = ( () -> entry.rule(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments...) )
        end
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
