export  calculateMessage!,
        calculateForwardMessage!,
        calculateBackwardMessage!,
        execute,
        clearMessages!

function execute(schedule_entry::ScheduleEntry)
    # Calculate the outbound message based on the inbound messages and the message calculation rule.
    # The resulting message is stored in the specified interface and is returned.

    outbound_interface = schedule_entry.node.interfaces[schedule_entry.outbound_interface_id]

    # Evaluate message calculation rule
    isdefined(schedule_entry, :execute) || error("Execute function not defined for schedule entry; maybe you forgot to prepare the schedule?")
    outbound_dist = schedule_entry.execute()

    # Post processing?
    if isdefined(schedule_entry, :post_processing)
        post_processed_output = schedule_entry.post_processing(outbound_dist)
        if (typeof(post_processed_output) <: ProbabilityDistribution) == false
            # Wrap the output in a DeltaDistribution before packing it in a Message
            post_processed_output = DeltaDistribution(post_processed_output)
        end
        outbound_interface.message = Message(post_processed_output) # Assign new post-processed message to interface
        outbound_dist = outbound_interface.message.payload # Reset outbound to post-processed result
    end

    # Print output for debugging
    if verbose
        node = schedule_entry.node
        interface = node.interfaces[schedule_entry.outbound_interface_id]
        interface_handle = (handle(interface)!="") ? "($(handle(interface)))" : ""
        println(replace("$(schedule_entry.rule) on $(typeof(node)) $(interface.node.id) interface $(schedule_entry.outbound_interface_id) $(interface_handle)", "ForneyLab.", ""))
        if isdefined(schedule_entry, :inbound_types) && isdefined(schedule_entry, :outbound_type)
            println(replace("$(schedule_entry.inbound_types) -> Message{$(schedule_entry.outbound_type)}", "ForneyLab.", ""))
        end
        if isdefined(schedule_entry, :post_processing)
            println(replace("Post processing: $(schedule_entry.post_processing)", "ForneyLab.", ""))
        end
        println("Result: $(outbound_dist)")
    end

    return outbound_dist
end

# Execute schedules
function execute(schedule::Schedule)
    # Execute a message passing schedule
    !isempty(schedule) || error("Cannot execute an empty schedule")

    # Print table header for execution log
    if verbose
        println("\n\nExecution log")
        println("--------------------------------------------")
    end

    for i=1:length(schedule)
        if verbose
            println("$(i).")
        end
        execute(schedule[i])
    end

    # Return the last message in the schedule
    entry = schedule[end] 
    return entry.node.interfaces[entry.outbound_interface_id].message
end

function clearMessages!(node::Node)
    # Clear all outbound messages on the interfaces of node
    for interface in node.interfaces
        clearMessage!(interface)
    end
end

function clearMessages!(edge::Edge)
    # Clear all messages on an edge.
    clearMessage!(edge.head.message)
    clearMessage!(edge.tail.message)
end

function inferOutboundType!(entry::ScheduleEntry, node::Node, outbound_interface_id::Int64, inbound_types::Vector, allowed_rules::Vector{Function})
    # Find all compatible calculation rules for the SumProduct algorithm
    outbound_types = []
    for update_function in allowed_rules
        available_rules = methods(update_function, [typeof(node); Type{Val{outbound_interface_id}}; inbound_types; Any])
        for rule in available_rules
            push!(outbound_types, rule.sig.types[end]) # Push the found outbound type on the stack
        end
    end

    # The outbound outbound_types should contain just one element (there should be just one available update rule)
    if length(outbound_types) == 0
        error("No calculation rule available for inbound types $(inbound_types) on node $(node)")
    elseif length(outbound_types) > 1
        error("Multiple outbound type possibilities ($(outbound_types)) for inbound types $(inbound_types) on node $(node)")
    elseif outbound_types[1] == Any
        # The computation rule produces Any, which indicates that the node is parametrized by its outbound type (e.g. a TerminalNode)
        (typeof(node).parameters[1] <: ProbabilityDistribution) || error("$(typeof(node)) $(node.id) must be parametrized by a ProbabilityDistribution")
        entry.outbound_type = typeof(node).parameters[1]
    elseif outbound_types[1] <: ProbabilityDistribution
        # There is only one possible outbound type and it is a probability distribution
        entry.outbound_type = outbound_types[1]
    else
        error("Unknown output of message calculation rule: $(outbound_types[1]) for node $(node)")
    end

    return entry
end

function inferOutboundType!(entry::ScheduleEntry, node::CompositeNode, outbound_interface_id::Int64, ::Vector, ::Vector{Function})
    # Infer outbound type of composite node
    outbound_interface = node.interfaces[outbound_interface_id]

    # Check if there is already a computation rule defined for this outbound interface
    if !haskey(node.computation_rules, outbound_interface)
        # Try to automatically generate a sum-product algorithm
        clearMessages!(node.internal_graph)
        internal_outbound_interface = node.interfaceid_to_terminalnode[outbound_interface_id].interfaces[1].partner
        composite_algo = SumProduct(internal_outbound_interface)
        node.computation_rules[outbound_interface] = composite_algo
    end

    # Fetch the algorithm that corresponds with outbound interface (either preset or just constructed)
    composite_algo = node.computation_rules[outbound_interface]
    (typeof(composite_algo) <: SumProduct) || error("Only sum product algorithms are currently supported on composite nodes")

    if !isdefined(composite_algo.schedule[end], :outbound_type) # Are the distribution types already inferred?
        inferDistributionTypes!(composite_algo) # Infer types on algorithm coupled to outbound interface of the composite node
    end

    entry.outbound_type = composite_algo.schedule[end].outbound_type # The last entry in the sum-product schedule holds the composite outbound type for this outbound interface

    return entry
end