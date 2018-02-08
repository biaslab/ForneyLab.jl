function writeInitializationBlock(schedule::Schedule, interface_to_msg_idx::Dict{Interface, Int})
    code = ""

    # Collect outbound types from schedule
    outbound_types = Dict()
    for entry in schedule
        outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    # Find breaker types from schedule outbound types
    breaker_types = Dict()
    for entry in schedule
        if entry.msg_update_rule <: ExpectationPropagationRule
            partner = ultimatePartner(entry.interface)
            breaker_types[partner] = outbound_types[partner]
        end
    end

    if !isempty(breaker_types) # Initialization block is only required when breakers are present
        code *= "function init$(name)()\n\n"

        # Write message (breaker) initialization code
        code *= "messages = Array{Message}($n_messages)\n"
        for (breaker_site, breaker_type) in breaker_types
            msg_idx = interface_to_msg_idx[breaker_site]
            breaker_str = replace(string(family(breaker_type)),"ForneyLab.", "") # Remove module prefixes
            code *= "messages[$(msg_idx)] = Message(vague($(breaker_str)))\n"
        end

        code *= "\nreturn messages\n\n"
        code *= "end\n\n"
    end

    return code
end

function writeMessagePassingBlock(schedule::Schedule, interface_to_msg_idx::Dict{Interface, Int})
    code = ""
    
    for schedule_entry in schedule
        # Collect inbounds and assign message id
        inbounds = collectInbounds(schedule_entry, schedule_entry.msg_update_rule, interface_to_msg_idx)

        # Apply update rule
        rule_id = schedule_entry.msg_update_rule
        rule_str = split(string(rule_id),'.')[end] # Remove module prefixes
        inbounds_str = join(inbounds, ", ")
        msg_idx = interface_to_msg_idx[schedule_entry.interface]
        code *= "messages[$msg_idx] = rule$(rule_str)($inbounds_str)\n"
    end

    return code
end

function writeMarginalsComputationBlock(schedule::MarginalSchedule, interface_to_msg_idx::Dict{Interface, Int})
    code = ""

    for schedule_entry in schedule
        if schedule_entry.marginal_update_rule == Void
            iface = schedule_entry.interfaces[1]
            code *= "marginals[:$(schedule_entry.target.id)] = messages[$(interface_to_msg_idx[iface])].dist\n"
        elseif schedule_entry.marginal_update_rule == Product
            iface1 = schedule_entry.interfaces[1]
            iface2 = schedule_entry.interfaces[2]
            code *= "marginals[:$(schedule_entry.target.id)] = messages[$(interface_to_msg_idx[iface1])].dist * messages[$(interface_to_msg_idx[iface2])].dist\n"
        else
            # Collect inbounds for marginal computation
            inbounds = collectInbounds(schedule_entry, interface_to_msg_idx)

            # Apply marginal update rule
            rule_id = schedule_entry.marginal_update_rule
            rule_str = split(string(rule_id),'.')[end] # Remove module prefixes
            inbounds_str = join(inbounds, ", ")
            code *= "marginals[:$(schedule_entry.target.id)] = rule$(rule_str)($inbounds_str)\n"
        end
    end

    return code
end

"""
Construct the inbound code that computes the marginal for `entry`.
Returns a vector with inbounds that correspond with required interfaces.
"""
function collectInbounds(entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(first(entry.target.edges))
    local_cluster_ids = localRecognitionFactorization(entry.target.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.target.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbounds
end

function messagePassingAlgorithm(schedule::Schedule, marginal_schedule::MarginalSchedule; file::String="", name::String="")
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    n_messages = length(schedule)

    # Assign message numbers to each interface in the schedule
    schedule = ForneyLab.flatten(schedule) # Inline all internal message passing
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)

    code = ""
    code *= writeInitializationBlock(schedule, interface_to_msg_idx)
    code *= "function step$(name)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}($n_messages))\n\n"
    code *= writeMessagePassingBlock(schedule, interface_to_msg_idx)
    code *= "\n"
    code *= writeMarginalsComputationBlock(marginal_schedule, interface_to_msg_idx)
    code *= "\nreturn marginals\n\n"
    code *= "end"

    # Write to file
    if !isempty(file)
        write(file, code)
    end

    return code
end

"""
Depending on the origin of the Clamp node message,
contruct the outbound message code.
"""
function messageString{T<:VariateType}(node::Clamp{T})
    var_type_str = split(string(T),'.')[end] # Remove module prefixes
    if node in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data array
        buffer, idx = ForneyLab.current_graph.placeholders[node]
        if idx > 0
            str = "Message($(var_type_str), PointMass, m=data[:$buffer][$idx])"
        else
            str = "Message($(var_type_str), PointMass, m=data[:$buffer])"
        end
    else
        # Insert constant
        str = "Message($(var_type_str), PointMass, m=$(node.value))"
    end

    return str
end

"""
Depending on the origin of the Clamp node message,
contruct the marginal code.
"""
function marginalString{T<:VariateType}(node::Clamp{T})
    var_type_str = split(string(T),'.')[end] # Remove module prefixes
    if node in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data array
        buffer, idx = ForneyLab.current_graph.placeholders[node]
        if idx > 0
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=data[:$buffer][$idx])"
        else
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=data[:$buffer])"
        end
    else
        # Insert constant
        str = "ProbabilityDistribution($(var_type_str), PointMass, m=$(node.value))"
    end

    return str
end