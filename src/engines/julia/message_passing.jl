function messagePassingAlgorithm(schedule::Schedule, marginal_schedule::MarginalSchedule=MarginalScheduleEntry[])
    # Assign message numbers to each interface in the schedule
    schedule = ForneyLab.flatten(schedule) # Inline all internal message passing
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)

    rf_dict = Dict()
    rf_dict[:schedule] = writeMessagePassingBlock(schedule, interface_to_msg_idx)
    writeInitializationBlock!(rf_dict, schedule, interface_to_msg_idx)
    rf_dict[:marginals] = writeMarginalsComputationBlock(marginal_schedule, interface_to_msg_idx)

    return rf_dict
end

function writeMessagePassingBlock(schedule::Schedule, interface_to_msg_idx::Dict{Interface, Int})
    schedule_vect = Vector{Dict}(undef, length(schedule))
    for (msg_idx, schedule_entry) in enumerate(schedule)
        schedule_entry_dict = schedule_vect[msg_idx]

        # Collect inbounds and assign message id
        schedule_entry_dict[:schedule_index] = msg_idx
        schedule_entry_dict[:message_update_rule] = schedule_entry.msg_update_rule
        schedule_entry_dict[:inbounds] = collectInbounds(schedule_entry, schedule_entry.msg_update_rule, interface_to_msg_idx)
    end

    return schedule_vect
end

"""
Depending on the origin of the Clamp node message, contruct the inbound stucture
"""
function messageDict(node::Clamp{T}) where T<:VariateType
    inbound = Dict(:variate_type => T,
                   :dist_or_msg  => Message)
    if node in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data array
        (buffer, idx) = ForneyLab.current_graph.placeholders[node]
        inbound[:buffer_id] = buffer
        if idx > 0
            inbound[:buffer_index] = idx
        end
    else
        # Insert constant
        inbound[:value] = node.value
    end

    return inbound
end

"""
Depending on the origin of the Clamp node message,
contruct the marginal code.
"""
function marginalDict(node::Clamp{T}) where T<:VariateType
    inbound = Dict(:variate_type => T,
                   :dist_or_msg  => ProbabilityDistribution)
    if node in keys(ForneyLab.current_graph.placeholders)
        # Distribution comes from data array
        (buffer, idx) = ForneyLab.current_graph.placeholders[node]
        inbound[:buffer_id] = buffer
        if idx > 0
            inbound[:buffer_index] = idx
        end
    else
        # Insert constant
        inbound[:value] = node.value
    end

    return inbound
end

function writeInitializationBlock!(rf_dict::Dict, schedule::Schedule, interface_to_msg_idx::Dict{Interface, Int}, n_messages::Int, name::String)
    # Collect outbound types from schedule
    outbound_types = Dict()
    for entry in schedule
        outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    # Find breaker types and dimensions
    update_clamp_flag = false # Flag that tracks whether the update of a clamped variable is required
    for entry in schedule
        partner = ultimatePartner(entry.interface)
        if (entry.msg_update_rule <: ExpectationPropagationRule)
            breaker_idx = interface_to_msg_idx[partner]
            breaker_dict = rf_dict[:schedule][breaker_idx]
            writeBreaker!(breaker_dict, outbound_types[partner], ()) # Univariate only
        elseif isa(entry.interface.node, Nonlinear) && (entry.interface == entry.interface.node.interfaces[2]) && (entry.interface.node.g_inv == nothing)
            # Set initialization in case of a nonlinear node without given inverse 
            iface = ultimatePartner(entry.interface.node.interfaces[2])
            breaker_idx = interface_to_msg_idx[iface]
            breaker_dict = rf_dict[:schedule][breaker_idx]
            writeBreaker!(breaker_dict, outbound_types[iface], entry.interface.node.dims)
        elseif !(partner == nothing) && isa(partner.node, Clamp)
            update_clamp_flag = true # Signifies the need for creating a custom `step!` function for optimizing clamped variables
            iface = entry.interface
            breaker_idx = interface_to_msg_idx[iface]
            breaker_dict = rf_dict[:schedule][breaker_idx]
            writeBreaker!(breaker_dict, outbound_types[iface], size(partner.node.value))
        end
    end

    rf_dict[:optimize] = update_clamp_flag

    return rf_dict
end

function writeBreaker!(breaker_dict::Dict, family::Type, dimensionality::Tuple)
    breaker_dict[:initialize] = true
    breaker_dict[:dimensionality] = dimensionality
    if family == Union{Gamma, Wishart} # Catch special case
        if dimensionality == ()
            breaker_dict[:family] = ForneyLab.Gamma
        else
            breaker_dict[:family] = ForneyLab.Wishart
        end
    else
        breaker_dict[:family] = family
    end

    return breaker_dict
end

function writeMarginalsComputationBlock(schedule::MarginalSchedule, interface_to_msg_idx::Dict{Interface, Int})
    schedule_vect = Vector{Dict}(undef, length(schedule))
    for (marg_idx, schedule_entry) in enumerate(schedule)
        if schedule_entry.marginal_update_rule == Nothing
            iface = schedule_entry.interfaces[1]
            inbounds = [Dict(:schedule_index => interface_to_msg_idx[iface])]
        elseif schedule_entry.marginal_update_rule == Product
            iface1 = schedule_entry.interfaces[1]
            iface2 = schedule_entry.interfaces[2]
            inbounds = [Dict(:schedule_index => interface_to_msg_idx[iface1]),
                        Dict(:schedule_index => interface_to_msg_idx[iface2])]
        else
            inbounds = collectInbounds(schedule_entry, interface_to_msg_idx)
        end

        schedule_vect[marg_idx] = Dict(:marginal_update_rule => schedule_entry.marginal_update_rule,
                                       :marginal_id => schedule_entry.target.id,
                                       :inbounds => inbounds)
    end

    return schedule_vect
end

"""
Construct the inbound code that computes the marginal for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
collectInbounds(entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int}) = collectMarginalNodeInbounds(entry.target.node, entry, interface_to_msg_idx)

function collectMarginalNodeInbounds(::FactorNode, entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = Dict[]
    entry_recognition_factor_id = recognitionFactorId(first(entry.target.edges))
    local_cluster_ids = localRecognitionFactorization(entry.target.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.target.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalDict(partner_node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, Dict(:schedule_index => inbound_idx))
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, Dict(:marginal_id => marginal_idx))
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbounds
end

