function messagePassingAlgorithm(schedule::Schedule)
    n_messages = length(schedule)

    code = "(data::Dict) -> begin\n\n"
    code *= "messages = Array{Message}($n_messages)\n"

    interface_to_msg_idx = Dict{Interface, Int}()
    for (msg_idx, schedule_entry) in enumerate(schedule)
        interface_to_msg_idx[schedule_entry.interface] = msg_idx
        node = schedule_entry.interface.node

        if isa(node, Constant)
            if node in keys(ForneyLab.current_graph.placeholders)
                # Read from data array
                buffer, idx = ForneyLab.current_graph.placeholders[node]
                if idx > 0
                    code *= "messages[$msg_idx] = Message{PointMass}(data[:$buffer][$idx])\n"
                else
                    code *= "messages[$msg_idx] = Message{PointMass}(data[:$buffer])\n"
                end
            else
                # Insert constant
                code *= "messages[$msg_idx] = Message{PointMass}($(node.value))\n"
            end
        else
            # Apply message update rule

            # Collect inbounds
            inbounds = String[]
            for node_interface in node.interfaces
                inbound_interface = node_interface.partner
                if haskey(interface_to_msg_idx, inbound_interface)
                    inbound_idx = interface_to_msg_idx[inbound_interface]
                    push!(inbounds, "messages[$inbound_idx]")
                else
                    push!(inbounds, "nothing")
                end
            end

            # Apply update rule
            rule_id = schedule_entry.msg_update_rule
            inbounds_str = join(inbounds, ", ")
            code *= "messages[$msg_idx] = rule($rule_id, [$inbounds_str])\n"
        end
    end

    # TODO: define ForneyLab.Julia.MessagePassing()
    # TODO: Constant and Placeholder are special cases for which the values
    # should be read from the data dictionary or node.value respectively

    code *= "\nend"

    return code
end
