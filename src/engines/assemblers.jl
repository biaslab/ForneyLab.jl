function assembleAlgorithm!(algo::Algorithm)
    # Generate structures for fast lookup
    algo.interface_to_schedule_entry = interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = targetToMarginalEntry(algo)

    for (id, rf) in algo.recognition_factorization
        rf.algorithm_id = algo.id
        assembleRecognitionFactor!(rf)
    end

    return algo
end

function assembleRecognitionFactor!(rf::RecognitionFactor)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    assembleMarginalTable!(rf)

    return rf
end

function assembleSchedule!(rf::RecognitionFactor)
    # Collect inbounds and assign message index per schedule entry
    for (msg_idx, schedule_entry) in enumerate(rf.schedule)
        schedule_entry.initialize = false # Preset initialization to false
        schedule_entry.inbounds = collectInbounds(schedule_entry, schedule_entry.message_update_rule)
        schedule_entry.schedule_index = msg_idx
    end

    return rf
end

function assembleInitialization!(rf::RecognitionFactor)
    interface_to_schedule_entry = current_algorithm.interface_to_schedule_entry

    # Collect outbound types from schedule
    outbound_types = Dict{Interface, Type}()
    for entry in rf.schedule
        outbound_types[entry.interface] = outboundType(entry.message_update_rule)
    end

    # Find breaker types and dimensions
    rf_update_clamp_flag = false # Flag that tracks whether the update of a clamped variable is required
    rf_initialize_flag = false # Indicate need for an initialization block
    for entry in rf.schedule
        partner = ultimatePartner(entry.interface)
        if (entry.message_update_rule <: ExpectationPropagationRule)
            breaker_entry = interface_to_schedule_entry[partner]
            assembleBreaker!(breaker_entry, family(outbound_types[partner]), ()) # Univariate only
            rf_initialize_flag = true 
        elseif isa(entry.interface.node, Nonlinear) && (entry.interface == entry.interface.node.interfaces[2]) && (entry.interface.node.g_inv == nothing)
            # Set initialization in case of a nonlinear node without given inverse 
            iface = ultimatePartner(entry.interface.node.interfaces[2])
            breaker_entry = interface_to_schedule_entry[iface]
            assembleBreaker!(breaker_entry, family(outbound_types[iface]), entry.interface.node.dims)
            rf_initialize_flag = true
        elseif !(partner == nothing) && isa(partner.node, Clamp)
            rf_update_clamp_flag = true # Signifies the need for creating a custom `step!` function for optimizing clamped variables
            iface = entry.interface
            breaker_entry = interface_to_schedule_entry[iface]
            assembleBreaker!(breaker_entry, family(outbound_types[iface]), size(partner.node.value))
            rf_initialize_flag = true
        end
    end

    rf.optimize = rf_update_clamp_flag
    rf.initialize = rf_initialize_flag

    return rf
end

function assembleMarginalTable!(rf::RecognitionFactor)
    interface_to_schedule_entry = current_algorithm.interface_to_schedule_entry
    
    for entry in rf.marginal_table
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

    return rf
end

function assembleBreaker!(breaker_entry::ScheduleEntry, family::Type, dimensionality::Tuple)
    breaker_entry.initialize = true
    breaker_entry.dimensionality = dimensionality
    if family == Union{Gamma, Wishart} # Catch special case
        if dimensionality == ()
            breaker_entry.family = ForneyLab.Gamma
        else
            breaker_entry.family = ForneyLab.Wishart
        end
    else
        breaker_entry.family = family
    end

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
        if idx > 0
            inbound.buffer_index = idx
        end
    end

    return inbound
end