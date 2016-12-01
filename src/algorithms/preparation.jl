###########################################
# Shared methods for algorithm construction
###########################################

function collectTypeVarValues(method::Method, argument_types::Vector{DataType}, results::Dict{Symbol,Any}=Dict{Symbol,Any}())
    # Collect the values of TypeVar objects in method.sig.types, based on the types of the actual arguments (argument_types)
    # The results are collected in a Dict, indexed by the names of the TypeVars.
    # Example: if method.sig.types[2] = Message{MvGaussian{dims}}
    #           and argument_types[1] = Message{MvGaussian{3}}
    #          {:dims => 3} is added to the results Dict.
    for a=5:length(method.sig.types) # skip the first 4 arguments since they do not contain inbounds (function, node, outbound_id, outbound_type)
        _collectTypeVarValues!(method.sig.types[a], argument_types[a-1], results) # Use internal function for recursive implementation
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
    if length(method_sig_type.parameters) != length(arg_type.parameters)
        error("Mismatch in parameter vector length between method: $(method_sig_type), and argument: $(arg_type). Consider implementing a custom _collectTypeVarValues! method")
    end

    if length(arg_type.parameters) > 0
        for p = 1:length(arg_type.parameters)
            _collectTypeVarValues!(method_sig_type.parameters[p], arg_type.parameters[p], results)
        end
    end
end

# Custom methods for matching distribution type parameters
function _collectTypeVarValues!(method_sig_type::Type{Univariate}, arg_type::Type{MvDelta{Float64}}, results::Dict{Symbol,Any})
    # Do nothing
end

function _collectTypeVarValues!{dims}(method_sig_type::Type{Multivariate{dims}}, arg_type::Type{MvDelta{Float64, dims}}, results::Dict{Symbol,Any})
    _collectTypeVarValues!(method_sig_type.parameters[1], arg_type.parameters[2], results)
end

function _collectTypeVarValues!{dims_n, dims_m}(method_sig_type::Type{MatrixVariate{dims_n, dims_m}}, arg_type::Type{MatrixDelta{Float64, dims_n, dims_m}}, results::Dict{Symbol,Any})
    _collectTypeVarValues!(method_sig_type.parameters[1], arg_type.parameters[2], results)
    _collectTypeVarValues!(method_sig_type.parameters[2], arg_type.parameters[3], results)
end

function _collectTypeVarValues!{dims}(method_sig_type::Type{MatrixVariate{dims, dims}}, arg_type::Type{Wishart{dims}}, results::Dict{Symbol,Any})
    _collectTypeVarValues!(method_sig_type.parameters[1], arg_type.parameters[1], results)
end

function collectTypeVarNames(dtype::DataType, results=Set{Symbol}())
    # Collect all names of TypeVar parameters in dtype
    for param in dtype.parameters
        if isa(param, TypeVar)
            push!(results, param.name)
        elseif isa(param, DataType)
            collectTypeVarNames(param, results) # recurse into the parameters of this parameter
        end
    end

    return results
end

function resolveTypeVars(dtype::TypeVar, method::Method, argument_types::Vector{DataType}, node::Node)
    # Return the value of dtype.
    # The value is determined from the method specification and the types of the actual arguments.
    typevar_values = collectTypeVarValues(method, argument_types)
    try
        return typevar_values[dtype.name]
    catch
        # Resort to node-specific fallback lookup.
        return outboundParameterValue(node, Val{typevar}, method, argument_types)
    end
end

function resolveTypeVars(dtype::DataType, method::Method, argument_types::Vector{DataType}, node::Node)
    # Build a DataType based on dtype, but with all TypeVar parameters replaced by their corresponding values.
    # The values are determined from the method specification and the types of the actual arguments.
    typevars = collectTypeVarNames(dtype)
    (length(typevars) > 0) || return dtype # Nothing to resolve
    typevar_values = collectTypeVarValues(method, argument_types) # Collect all TypeVar values in a Dict indexed by the TypeVar name
    for typevar in typevars
        if !haskey(typevar_values, typevar)
            # Call node-specific function to resolve typevar
            typevar_values[typevar] = outboundParameterValue(node, Val{typevar}, method, argument_types)
        end
    end

    return _resolveTypeVars(dtype, typevar_values) # Build the resolved DataType
end

function _resolveTypeVars(dtype::DataType, lookup_table::Dict{Symbol,Any})
    # Internal function called by resolveTypeVars.
    # Returns dtype with all TypeVar parameters replaced by the values defined in lookup_table.

    if length(dtype.parameters) == 0
        return dtype
    else
        param_values = map(dt -> _resolveTypeVars(dt, lookup_table), dtype.parameters)
        return Main.eval(parse("$(dtype.name){" * join(param_values,",") * "}"))
    end
end

_resolveTypeVars(dtype::TypeVar, lookup_table::Dict{Symbol,Any}) = lookup_table[dtype.name]

function collectAllOutboundTypes(rule::Function, call_signature::Vector, node::Node)
    # Collect all outbound distribution types that can be generated with the specified rule and calling_signature combination.
    # The node argument is only meant for overloading purposes (see TerminalNode for example).
    # Returns: Vector{DataType} containing the possible outbound types.
    # The entries of the returned vector can be <: ProbabilityDistribution or Approximation

    outbound_types = DataType[]

    for method in methods(rule, call_signature)
        ob_type = resolveTypeVars(method.sig.types[4], method, call_signature, node)

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
            return DataType[Partitioned{ob_type, num_factors} for ob_type in collectAllOutboundTypes(rule, sig, node)]
        else
            return DataType[]
        end
    else
        sig = buildRuleSignature(rule, node, outbound_interface_id, inbound_types, approx)
        return collectAllOutboundTypes(rule, sig, node)
    end
end

function inferOutboundType!{rule}(entry::ScheduleEntry{rule})
    # Infer the outbound type from the node type and the types of the inbounds

    inbound_types = entry.inbound_types
    rule_implementation = implementation(rule)
    if isdefined(entry, :outbound_type)
        # The outbound type is already fixed, so we just need to validate that there exists a suitable computation rule

        # Try to find a matching exact rule
        for outbound_type in collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types)
            (entry.outbound_type <: outbound_type) && return(entry)
        end

        # Try to find a matching exact rule by unrolling partitioned inbounds
        for outbound_type in collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, unroll_partitioned_inbounds=true)
            if entry.outbound_type <: outbound_type
                entry.unrolling_factor = outbound_type.parameters[2]
                return(entry)
            end
        end

        # Try to find a matching approximate rule
        approx = isdefined(entry, :approximation) ? Type{entry.approximation} : Any
        for outbound_type in collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, approx)
            if entry.outbound_type <: outbound_type.parameters[1]
                entry.approximation = outbound_type.parameters[2]
                return entry
            end
        end

        # Try to find a matching approximate rule by unrolling partitioned inbounds
        for outbound_type in collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, approx, unroll_partitioned_inbounds=true)
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
        outbound_types = collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type.")
        end

        # Try to apply an exact rule by unrolling partitioned inbounds
        outbound_types = collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, unroll_partitioned_inbounds=true)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1]
            entry.unrolling_factor = outbound_types[1].parameters[2]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type.")
        end

        # Try to apply an approximate rule
        approx = isdefined(entry, :approximation) ? Type{entry.approximation} : Any
        outbound_types = collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, approx)
        if length(outbound_types) == 1
            entry.outbound_type = outbound_types[1].parameters[1]
            entry.approximation = outbound_types[1].parameters[2]
            return entry
        elseif length(outbound_types) > 1
            error("There are multiple outbound type possibilities for schedule entry:\n$(entry)Inbound types: $(inbound_types)\nPlease specify a message type and if required also an approximation method.")
        end

        # Try to apply an approximate rule by unrolling partitioned inbounds
        outbound_types = collectOutboundTypes(rule_implementation, entry.node, entry.outbound_interface_id, inbound_types, approx, unroll_partitioned_inbounds=true)
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
        if ib_type <: Partitioned
            return ib_type.parameters[2]
        elseif ((ib_type <: Message) && (ib_type.parameters[1] <: Partitioned))
            return ib_type.parameters[1].parameters[2]
        end
    end

    return 0
end

function stripPartitionedInboundTypes(inbound_types::Vector{DataType})
    # Simplify a list of inbound types that contains Partitioned or Message{Partitioned} entries.
    # Replace every partitioned inbound type by the type of the factors.
    # Check if all partitioned inbounds contain the same number of factors. If not, throw an error.
    stripped_inbound_types = Vector{DataType}()
    num_factors = 0
    for ib_type in inbound_types
        if ib_type <: Partitioned
            part_type = ib_type
        elseif ((ib_type <: Message) && (ib_type.parameters[1] <: Partitioned))
            part_type = ib_type.parameters[1]
        else
            # Regular inbound type
            push!(stripped_inbound_types, ib_type)
            continue
        end

        if num_factors == 0
            num_factors = part_type.parameters[2]
        elseif num_factors != part_type.parameters[2]
            error("Partitioned inbounds for a message computation rule should have the same number of factors.")
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

function buildExecute!{rule}(   entry::ScheduleEntry{rule},
                                inbound_arguments::Vector;
                                outbound_dist::ProbabilityDistribution=entry.node.interfaces[entry.outbound_interface_id].message.payload)
    # Construct the entry.execute function.
    # This function is called by the prepare methods of inference algorithms.

    # Get pointer to the outbound distribution
    rule_implementation = implementation(rule)

    # Save the "compiled" message computation rule as an anomynous function in entry.execute
    if entry.unrolling_factor > 0
        # Unroll partitioned inbounds and call message calculation rule for each factor
        factor_inbounds = copy(inbound_arguments) # inbounds for a single factor
        inbounds_to_unroll = Int64[]
        for i=1:length(inbound_arguments)
            if typeof(inbound_arguments[i]) <: Partitioned
                factor_inbounds[i] = inbound_arguments[i].factors[1]
                push!(inbounds_to_unroll, i)
            elseif ((typeof(inbound_arguments[i]) <: Message) && (typeof(inbound_arguments[i].payload) <: Partitioned))
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
                    rule_implementation(entry.node, Val{entry.outbound_interface_id}, outbound_dist.factors[factor_idx], factor_inbounds..., entry.approximation)
                else
                    rule_implementation(entry.node, Val{entry.outbound_interface_id}, outbound_dist.factors[factor_idx], factor_inbounds...)
                end
            end

            return outbound_dist
        end
    else
        if isdefined(entry, :approximation)
            entry.execute = ( () -> rule_implementation(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments..., entry.approximation) )
        else
            entry.execute = ( () -> rule_implementation(entry.node, Val{entry.outbound_interface_id}, outbound_dist, inbound_arguments...) )
        end
    end

    return entry
end

function injectParameters!{T<:ProbabilityDistribution}(destination::T, source::T)
    # Fill the parameters of a destination distribution with the copied parameters of the source

    for field in fieldnames(source)
        if isdefined(source, field)
            setfield!(destination, field, deepcopy(getfield(source, field)))
        end
    end

    return destination
end
