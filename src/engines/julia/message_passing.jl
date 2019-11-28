function assembleAlgorithm!(rf::RecognitionFactor)
    # Assign message numbers to each interface in the schedule
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    target_to_marginal_entry = ForneyLab.targetToMarginalEntry(rf.marginal_schedule)

    assembleSchedule!(rf.schedule, interface_to_schedule_entry, target_to_marginal_entry)
    assembleInitialization!(rf, interface_to_schedule_entry)
    assembleMarginalSchedule!(rf.marginal_schedule, interface_to_schedule_entry, target_to_marginal_entry)

    return rf
end

function assembleSchedule!(schedule::Schedule, interface_to_schedule_entry::Dict, target_to_marginal_entry::Dict)
    # Collect inbounds and assign message index per schedule entry
    for (msg_idx, schedule_entry) in enumerate(schedule)
        schedule_entry.inbounds = collectInbounds(schedule_entry, schedule_entry.message_update_rule, interface_to_schedule_entry, target_to_marginal_entry)
        schedule_entry.schedule_index = msg_idx
    end

    return schedule
end

function assembleInitialization!(rf::RecognitionFactor, interface_to_schedule_entry::Dict)
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

function assembleMarginalSchedule!(schedule::MarginalSchedule, interface_to_schedule_entry::Dict, target_to_marginal_entry::Dict)
    for entry in schedule
        if entry.marginal_update_rule == Nothing
            iface = entry.interfaces[1]
            inbounds = [interface_to_schedule_entry[iface]]
        elseif entry.marginal_update_rule == Product
            iface1 = entry.interfaces[1]
            iface2 = entry.interfaces[2]
            inbounds = [interface_to_schedule_entry[iface1], 
                        interface_to_schedule_entry[iface2]]
        else
            inbounds = collectInbounds(entry, interface_to_schedule_entry, target_to_marginal_entry)
        end

        entry.marginal_id = entry.target.id
        entry.inbounds = inbounds
    end

    return schedule
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

"""
Construct the inbound code that computes the marginal for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
collectInbounds(entry::MarginalScheduleEntry, interface_to_schedule_entry::Dict, target_to_marginal_entry::Dict) = collectMarginalNodeInbounds(entry.target.node, entry, interface_to_schedule_entry, target_to_marginal_entry)

function collectMarginalNodeInbounds(::FactorNode, entry::MarginalScheduleEntry, interface_to_schedule_entry::Dict, target_to_marginal_entry::Dict)
    # Collect inbounds
    inbounds = Any[]
    entry_recognition_factor = recognitionFactor(first(entry.target.edges))
    local_clusters = localRecognitionFactorization(entry.target.node)

    recognition_factors = Union{RecognitionFactor, Edge}[] # Keep track of encountered recognition factors
    for node_interface in entry.target.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor = recognitionFactor(node_interface.edge)

        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(partner_node, ProbabilityDistribution))
        elseif node_interface_recognition_factor === entry_recognition_factor
            # Collect message from previous result
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif !(node_interface_recognition_factor in recognition_factors)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            target = local_clusters[node_interface_recognition_factor]
            push!(inbounds, target_to_marginal_entry[target])
        end

        push!(recognition_factors, node_interface_recognition_factor)
    end

    return inbounds
end

