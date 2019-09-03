function writeInitializationBlock(schedule::Schedule, interface_to_msg_idx::Dict{Interface, Int}, n_messages::Int, name::String)
    code = ""

    # Collect outbound types from schedule
    outbound_types = Dict()
    for entry in schedule
        outbound_types[entry.interface] = outboundType(entry.msg_update_rule)
    end

    # Find breaker types and dimensions
    update_clamp_flag = false # Flag that tracks whether the update of a clamped variable is required
    breaker_types = Dict()
    breaker_dims = Dict()
    for entry in schedule
        partner = ultimatePartner(entry.interface)
        if (entry.msg_update_rule <: ExpectationPropagationRule)
            breaker_types[partner] = outbound_types[partner]
            breaker_dims[partner] = () # Univariate only
        elseif isa(entry.interface.node, Nonlinear) && (entry.interface == entry.interface.node.interfaces[2]) && (entry.interface.node.g_inv == nothing)
            # Set initialization in case of a nonlinear node without given inverse 
            iface = ultimatePartner(entry.interface.node.interfaces[2])
            breaker_types[iface] = outbound_types[iface]
            breaker_dims[iface] = entry.interface.node.dims
        elseif !(partner == nothing) && isa(partner.node, Clamp)
            update_clamp_flag = true # Signifies the need for creating a custom `step!` function for optimizing clamped variables
            iface = entry.interface
            breaker_types[iface] = outbound_types[iface]
            breaker_dims[iface] = size(partner.node.value)
        end
    end

    if update_clamp_flag # Explain the need for a custom `step!` definition
        code *= "# You have created an algorithm that requires updates for (a) clamped parameter(s).\n"
        code *= "# This algorithm requires the definition of a custom `optimize!` function that updates the parameter value(s)\n"
        code *= "# by altering the `data` dictionary in-place. The custom `optimize!` function may be based on the mockup below:\n\n"
        code *= "# function optimize$(name)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init$(name)())\n"
        code *= "# \t...\n"
        code *= "# \treturn data\n"
        code *= "# end\n\n"
    end

    if !isempty(breaker_types) # Initialization block is only required when breakers are present
        code *= "function init$(name)()\n\n"

        # Write message (breaker) initialization code
        code *= "messages = Array{Message}(undef, $n_messages)\n"
        for (breaker_site, breaker_type) in breaker_types
            msg_idx = interface_to_msg_idx[breaker_site]
            breaker_type_str = replace(string(family(breaker_type)),"ForneyLab." => "") # Remove module prefixes
            if breaker_dims[breaker_site] == () # Univariate
                code *= "messages[$(msg_idx)] = Message(vague($(breaker_type_str)))\n"
            else
                code *= "messages[$(msg_idx)] = Message(vague($(breaker_type_str), $(breaker_dims[breaker_site])))\n"
            end
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
        if schedule_entry.marginal_update_rule == Nothing
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

    !isempty(schedule) && (code *= "\n")

    return code
end

"""
Construct the inbound code that computes the marginal for `entry`. Allows for
overloading and for a user the define custom node-specific inbounds collection.
Returns a vector with inbounds that correspond with required interfaces.
"""
collectInbounds(entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int}) = collectMarginalNodeInbounds(entry.target.node, entry, interface_to_msg_idx)

function collectMarginalNodeInbounds(::FactorNode, entry::MarginalScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
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

function messagePassingAlgorithm(schedule::Schedule, marginal_schedule::MarginalSchedule=MarginalScheduleEntry[]; file::String="", name::String="")
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    n_messages = length(schedule)

    # Assign message numbers to each interface in the schedule
    schedule = ForneyLab.flatten(schedule) # Inline all internal message passing
    schedule = ForneyLab.condense(schedule) # Remove Clamp node entries
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)

    code = ""
    code *= writeInitializationBlock(schedule, interface_to_msg_idx, n_messages, name)
    code *= "function step$(name)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, $n_messages))\n\n"
    code *= writeMessagePassingBlock(schedule, interface_to_msg_idx)
    code *= "\n"
    code *= writeMarginalsComputationBlock(marginal_schedule, interface_to_msg_idx)
    code *= "return marginals\n\n"
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
function messageString(node::Clamp{T}) where T<:VariateType
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
        if size(node.value) == (1,1)
            str = "Message($(var_type_str), PointMass, m=mat($(node.value[1])))"
        elseif isa(node.value, Diagonal)
            str = "Message($(var_type_str), PointMass, m=Diagonal($(node.value.diag)))"
        elseif isa(node.value, AbstractMatrix) && (size(node.value)[2] == 1)
            str = "Message($(var_type_str), PointMass, m=[$(join(node.value, " "))]')"
        else
            str = "Message($(var_type_str), PointMass, m=$(node.value))"
        end
    end

    return str
end

"""
Depending on the origin of the Clamp node message,
contruct the marginal code.
"""
function marginalString(node::Clamp{T}) where T<:VariateType
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
        if size(node.value) == (1,1)
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=mat($(node.value[1])))"
        elseif isa(node.value, Diagonal)
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=Diagonal($(node.value.diag)))"
        elseif isa(node.value, AbstractMatrix) && (size(node.value)[2] == 1)
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=[$(join(node.value, " "))]')"
        else
            str = "ProbabilityDistribution($(var_type_str), PointMass, m=$(node.value))"
        end
    end

    return str
end