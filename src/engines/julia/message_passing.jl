export load, freeEnergyAlgorithm

function freeEnergyAlgorithm(recognition_factors::Vector{RecognitionFactor}=collect(values(current_recognition_factorization.recognition_factors)))
    # Collect nodes connected to external edges
    nodes_connected_to_external_edges = Set{FactorNode}()
    for rf in recognition_factors
        union!(nodes_connected_to_external_edges, nodesConnectedToExternalEdges(rf))
    end

    # Write evaluation function for free energy
    energy_block = ""
    entropy_block = ""
    for node in sort(collect(nodes_connected_to_external_edges))
        # Average energy
        inbounds = collectMarginals(node)
        inbounds_str = join(inbounds, ", ")
        node_str = replace(string(typeof(node)),"ForneyLab.", "") # Remove module prefixes
        energy_block *= "F += averageEnergy($(node_str), $(inbounds_str))\n"

        # Differential entropy
        partner_iface = ultimatePartner(node.interfaces[1])
        if !(partner_iface == nothing) && !isa(partner_iface.node, Clamp)
            entropy_block *= "F -= differentialEntropy(marginals[:$(partner_iface.edge.variable.id)])\n"
        end
    end

    # Combine blocks
    code = "function freeEnergy(data::Dict, marginals::Dict)\n\n"
    code *= "F = 0.0\n\n"
    code *= energy_block*"\n"entropy_block
    code *= "\nreturn F\n"
    code *= "end"

    return code
end

"""
Collect marginals associated with all edges connected to `node`
"""
function collectMarginals(node::FactorNode)
    # Collect marginals
    inbound_marginals = String[]
    for interface in node.interfaces
        partner_node = ultimatePartner(interface).node
        if isa(partner_node, Clamp)
            # Hard-code marginal of constant node
            push!(inbound_marginals, marginalString(partner_node))
        else
            # Collect marginal from marginal dictionary
            push!(inbound_marginals, "marginals[:$(interface.edge.variable.id)]")
        end
    end

    return inbound_marginals
end

# TODO: in-place operations for message and marginal computations?
function messagePassingAlgorithm(schedule::Schedule, targets::Vector{Variable}=Variable[]; file::String="", name::String="")
    schedule = ForneyLab.flatten(schedule) # Inline all internal message passing
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    n_messages = length(schedule)

    # Assign message numbers to each interface in the schedule
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)

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
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
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
collectInbounds{T<:VariationalRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int}) = collectVariationalNodeInbounds(entry.interface.node, entry, interface_to_msg_idx)

function collectVariationalNodeInbounds(::FactorNode, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge)
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        partner_node = inbound_interface.node
        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, "nothing")
        elseif isa(partner_node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(partner_node))
        elseif recognitionFactorId(node_interface.edge) == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        else
            # Collect marginal from marginal dictionary
            push!(inbounds, "marginals[:$(node_interface.edge.variable.id)]")
        end
    end

    return inbounds
end

"""
Collect and construct EP update code for each inbound.
"""
function collectInbounds{T<:ExpectationPropagationRule}(entry::ScheduleEntry, ::Type{T}, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
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