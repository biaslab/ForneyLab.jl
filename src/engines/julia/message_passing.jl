export load

# TODO: in-place operations for message and marginal computations?
function messagePassingAlgorithm(schedule::Schedule, targets::Vector{Variable}=Variable[]; file::String="", name::String="")
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    n_messages = length(schedule)

    # Assign message numbers to each interface in the schedule
    interface_to_msg_idx = Dict{Interface, Int}()
    for (msg_idx, schedule_entry) in enumerate(schedule)
        interface_to_msg_idx[schedule_entry.interface] = msg_idx
    end

    # Collect outbound types from schedule
    outbound_types = Dict()
    for entry in schedule
        outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    code = ""


    #--------------------------------
    # Write initialization code block
    #--------------------------------

    # Find breaker types from schedule outbound types
    breaker_types = Dict()
    for entry in schedule
        if entry.msg_update_rule <: ExpectationPropagationRule
            breaker_types[entry.interface.partner] = outbound_types[entry.interface.partner]
        end
    end

    if !isempty(breaker_types) # Initialization block is only required when breakers are present
        code *= "function init$(name)()\n\n"
        
        # Write message (breaker) initialization code
        code *= "messages = Array{Message}($n_messages)\n"
        for (breaker_site, breaker_type) in breaker_types
            msg_idx = interface_to_msg_idx[breaker_site]
            breaker_str = replace(string(breaker_type),"ForneyLab.", "") # Remove module prefixes
            code *= "messages[$(msg_idx)] = vague($(breaker_str))\n"
        end

        code *= "\nreturn messages\n\n"
        code *= "end\n\n"
    end


    #---------------------------------
    # Write message passing code block
    #---------------------------------

    code *= "function step$(name)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}($n_messages))\n\n"

    # Write message computation code
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

    # Write marginal computation code
    code *= "\n"
    for variable in targets
        target_edge = first(variable.edges) # For the sake of consistency, we always take the first edge.
        if target_edge.a == nothing # Handle cases where there is a `dangling` edge
            msg_id_b = interface_to_msg_idx[target_edge.b]
            code *= "marginals[:$(variable.id)] = messages[$msg_id_b].dist\n"
        elseif target_edge.b == nothing
            msg_id_a = interface_to_msg_idx[target_edge.a]
            code *= "marginals[:$(variable.id)] = messages[$msg_id_a].dist\n"
        else
            msg_id_a = interface_to_msg_idx[target_edge.a]
            msg_id_b = interface_to_msg_idx[target_edge.b]
            code *= "marginals[:$(variable.id)] = messages[$msg_id_a].dist * messages[$msg_id_b].dist\n"
        end
    end

    code *= "\nreturn marginals\n\n"
    code *= "end"

    # Write to file
    if !isempty(file)
        write(file, code)
    end

    return code
end

messagePassingAlgorithm(schedule::Schedule, target::Variable; file::String="", name::String="") = messagePassingAlgorithm(schedule, [target], file=file, name=name)

"""
Collect and construct SP update code for each inbound.
"""
function collectInbounds{T<:SumProductRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    node = entry.interface.node
    for node_interface in node.interfaces
        inbound_interface = node_interface.partner
        if node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbound_messages, "nothing")
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    return inbound_messages
end

"""
Collect and construct VMP update code for each inbound.
"""
function collectInbounds{T<:VariationalRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbound_marginals = String[]
    node = entry.interface.node
    for node_interface in node.interfaces
        inbound_interface = node_interface.partner
        partner_node = inbound_interface.node
        if node_interface == entry.interface
            push!(inbound_marginals, "nothing") # Could also just push the marginal on the edge instead of "nothing"
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbound_marginals, marginalString(partner_node))
        else
            # Collect marginal from marginal dictionary
            push!(inbound_marginals, "marginals[:$(node_interface.edge.variable.id)]")
        end
    end

    return inbound_marginals
end

"""
Collect and construct EP update code for each inbound.
"""
function collectInbounds{T<:ExpectationPropagationRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    node = entry.interface.node
    for node_interface in node.interfaces
        inbound_interface = node_interface.partner
        if isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    return inbound_messages
end

"""
Depending on the origin of the Clamp node message,
contruct the outbound message code.
"""
function messageString(node::Clamp)
    if node in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data array
        buffer, idx = ForneyLab.current_graph.placeholders[node]
        if idx > 0
            str = "Message(PointMass, m=data[:$buffer][$idx])"
        else
            str = "Message(PointMass, m=data[:$buffer])"
        end
    else
        # Insert constant
        str = "Message(PointMass, m=$(node.value))"
    end

    return str
end

"""
Depending on the origin of the Clamp node message,
contruct the marginal code.
"""
function marginalString(node::Clamp)
    if node in keys(ForneyLab.current_graph.placeholders)
        # Message comes from data array
        buffer, idx = ForneyLab.current_graph.placeholders[node]
        if idx > 0
            str = "ProbabilityDistribution(PointMass, m=data[:$buffer][$idx])"
        else
            str = "ProbabilityDistribution(PointMass, m=data[:$buffer])"
        end
    else
        # Insert constant
        str = "ProbabilityDistribution(PointMass, m=$(node.value))"
    end

    return str
end