function assembleInferenceAlgorithm!(algo::InferenceAlgorithm)
    # Generate structures for fast lookup
    algo.interface_to_schedule_entry = interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = targetToMarginalEntry(algo)

    for (id, pf) in algo.posterior_factorization
        pf.algorithm_id = algo.id
        assemblePosteriorFactor!(pf)
    end

    return algo
end

function assemblePosteriorFactor!(pf::PosteriorFactor)
    assembleSchedule!(pf)
    assembleInitialization!(pf)
    assembleMarginalTable!(pf)

    return pf
end

function assembleSchedule!(pf::PosteriorFactor)
    # Collect inbounds and assign message index per schedule entry
    for (msg_idx, schedule_entry) in enumerate(pf.schedule)
        schedule_entry.initialize = false # Preset initialization to false
        schedule_entry.inbounds = collectInbounds(schedule_entry, schedule_entry.message_update_rule)
        schedule_entry.schedule_index = msg_idx
    end

    return pf
end

function assembleInitialization!(pf::PosteriorFactor)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry

    # Collect outbound types from schedule
    outbound_types = Dict{Interface, Type}()
    for entry in pf.schedule
        outbound_types[entry.interface] = outboundType(entry.message_update_rule)
    end

    # Find breaker types and dimensions
    pf_update_clamp_flag = false # Flag that tracks whether the update of a clamped variable is required
    pf_initialize_flag = false # Indicate need for an initialization block
    for entry in pf.schedule
        partner = ultimatePartner(entry.interface)

        # Assemble breakers
        if entry.interface in pf.breaker_interfaces
            (_, dims) = breakerParameters(entry.interface)
            assembleBreaker!(entry, family(outbound_types[entry.interface]), dims)
            pf_initialize_flag = true # Signifies the need for an initialization block
        elseif !(partner == nothing) && isa(partner.node, Clamp)
            assembleBreaker!(entry, family(outbound_types[entry.interface]), size(partner.node.value))
            pf_initialize_flag = true
            pf_update_clamp_flag = true # Signifies the need for creating a custom `step!` function for optimizing clamped variables
        end
    end

    pf.optimize = pf_update_clamp_flag
    pf.initialize = pf_initialize_flag

    return pf
end

function assembleMarginalTable!(pf::PosteriorFactor)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    
    for entry in pf.marginal_table
        if entry.marginal_update_rule == Nothing
            iface = entry.interfaces[1]
            inbounds = [interface_to_schedule_entry[iface]]
        elseif entry.marginal_update_rule == Product
            iface1 = entry.interfaces[1]
            iface2 = entry.interfaces[2]
            inbounds = [interface_to_schedule_entry[iface1], 
                        interface_to_schedule_entry[iface2]]
        else
            inbounds = collectInbounds(entry)
        end

        entry.marginal_id = entry.target.id
        entry.inbounds = inbounds
    end

    return pf
end

function assembleBreaker!(breaker_entry::ScheduleEntry, family::Type, dimensionality::Tuple)
    breaker_entry.initialize = true
    breaker_entry.family = family
    breaker_entry.dimensionality = dimensionality

    return breaker_entry
end

"""
Depending on the origin of the Clamp node message, contruct a message or marginal inbound
"""
function assembleClamp!(inbound::Clamp, dist_or_msg::Type)
    inbound.dist_or_msg = dist_or_msg
    if inbound in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data buffer
        (buffer, idx) = ForneyLab.current_graph.placeholders[inbound]
        inbound.buffer_id = buffer
        inbound.buffer_index = idx # Can also be 0, in which case index is ignored
    end

    return inbound
end