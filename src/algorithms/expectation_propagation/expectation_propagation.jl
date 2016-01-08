export ExpectationPropagation

type ExpectationPropagation <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
    sites::Vector{Interface}
    num_iterations::Int64
    callback::Function
end

function ExpectationPropagation(
            sites::Vector{Interface};
            num_iterations::Int64 = 100,
            callback::Function = () -> false
            post_processing_functions = Dict{Interface, Function}())
    # Build an EP message passing algorithm for the specified sites.
    # num_iterations specifies the maximum number of iterations.
    # After each iteration, callback is called to allow for convergence checks.
    # If the callback returns true, the algorithm is terminated.

    # InferenceAlgorithm overview:
    # 1. Init all sites with vague messages
    # 2. For all sites i=1:N
    #   2a. Calculate cavity distribution i
    #   2b. Calculate site distribution i
    # 3. Check stopping criteria, goto 2

    (length(sites) > 0) || error("Specify at least one site")

    # Build schedule
    total_schedule = Vector{Interface}()
    for i = 1:length(sites)
        site = sites[i]
        # Prepend sites b/c of vague initialization
        total_schedule = vcat(sites, total_schedule)
        # Add schedule for cavity distribution to total_schedule
        total_schedule = generateScheduleByDFS!(site.partner, total_schedule)
        total_schedule = total_schedule[length(sites)+1:end] # Strip sites prepend
        # Build list of other sites, prepend to total schedule
        if i < length(sites)
            other_sites = vcat(sites[1:i-1], sites[i+1:end])
        else
            other_sites = sites[1:i-1]
        end
        total_schedule = vcat(other_sites, total_schedule)
        total_schedule = generateScheduleByDFS!(site, total_schedule)
        total_schedule = total_schedule[length(sites):end] # Strip other sites prepend
    end

    schedule = convert(Schedule, total_schedule, sumProduct!)
    for schedule_entry in schedule
        if schedule_entry.interface in sites
            schedule_entry.rule = ep!
        end
    end
    setPostProcessing!(schedule, post_processing_functions)

    # Build execute function
    function exec(algorithm)
        # Init all sites with vague messages
        for site in algorithm.sites
            vague!(site.message.payload)
        end
        # Execute schedule until stopping criterium is met
        for iteration_count = 1:algorithm.num_iterations
            execute(algorithm.schedule)
            # Check stopping criteria
            if algorithm.callback()
                break
            end
        end
    end

    algo = ExpectationPropagation(exec, schedule, sites, num_iterations, callback)
    inferDistributionTypes!(algo)

    return algo
end

############################################
# Type inference and preparation
############################################

function inferDistributionTypes!(algo::ExpectationPropagation)
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}() # Lookup table from interface to schedule entry

    for entry in algo.schedule
        collectInboundTypes!(entry, schedule_entries, algo) # Fill entry.inbound_types
        inferOutboundType!(entry, [sumProduct!, ep!]) # Infer the outbound message type, fill entry.outbound_type

        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        schedule_entries[outbound_interface] = entry # Add entry to lookup table
    end

    return algo
end

function collectInboundTypes!(entry::ScheduleEntry, schedule_entries::Dict{Interface, ScheduleEntry}, ::ExpectationPropagation)
    # Look up the types of the inbound messages for entry.
    # Fill entry.inbound_types
    entry.inbound_types = []

    for (id, interface) in enumerate(entry.node.interfaces)
        if (id == entry.outbound_interface_id) && (entry.rule == sumProduct!)
            # Incoming msg on outbound interface is always Void for sumProduct! rule
            push!(entry.inbound_types, Void)
        else
            push!(entry.inbound_types, Message{schedule_entries[interface.partner].outbound_type})
        end
    end

    return entry
end

function prepare!(algo::ExpectationPropagation)
    # Populate the graph with vague messages of the correct types
    for entry in algo.schedule
        ensureMessage!(entry.node.interfaces[entry.outbound_interface_id], entry.outbound_type)
    end

    # Compile the schedule
    compile!(algo.schedule, algo)

    return algo
end

function compile!(entry::ScheduleEntry, ::Type{Val{symbol(ep!)}}, ::InferenceAlgorithm)
    # Generate entry.execute for schedule entry with ep! calculation rule

    inbound_messages = [interface.partner.message for interface in entry.node.interfaces]

    return buildExecute!(entry, inbound_messages)
end