import Base.show
export ExpectationPropagation

type ExpectationPropagation <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
    sites::Vector{Interface}
    n_iterations::Int64
    callback::Function
end

function show(algo::ExpectationPropagation)
    println("ExpectationPropagation inference algorithm")
    println("    # sites: $(length(algo.sites))")
    println("    max. number of iterations: $(algo.n_iterations)")
    println("    callback function: $(algo.callback)")
    println("Use show(algo.schedule) to view the message passing schedule.")
end

function ExpectationPropagation(
            sites::Vector{Tuple{Interface, DataType}};
            n_iterations::Int64 = 100,
            callback::Function = ( () -> false ),
            post_processing_functions = Dict{Interface, Function}())

    # Build an EP message passing algorithm for the specified sites.
    # sites is a list of (interface, recognition_distribution) tuples,
    # where recognition_distribution <: ProbabilityDistribution.
    # n_iterations specifies the maximum number of iterations.
    # After each iteration, callback is called to allow for convergence checks.
    # If the callback returns true, the algorithm is terminated.

    # InferenceAlgorithm overview:
    # 1. Init all sites with vague messages
    # 2. For all sites i=1:N
    #   2a. Calculate cavity distribution i
    #   2b. Calculate site distribution i
    # 3. Check stopping criteria, goto 2

    (length(sites) > 0) || error("Specify at least one site")

    sitelist = Interface[site[1] for site in sites]
    recognition_distributions = Dict{Interface,DataType}([interface => distribution_type for (interface, distribution_type) in sites])

    # Build schedule
    total_schedule = Vector{Interface}()
    for i = 1:length(sitelist)
        site = sitelist[i]
        # Prepend sitelist b/c of vague initialization
        total_schedule = vcat(sitelist, total_schedule)
        # Add schedule for cavity distribution to total_schedule
        total_schedule = generateScheduleByDFS!(site.partner, total_schedule)
        total_schedule = total_schedule[length(sitelist)+1:end] # Strip sitelist prepend
        # Build list of other sites, prepend to total schedule
        if i < length(sitelist)
            other_sites = vcat(sitelist[1:i-1], sitelist[i+1:end])
        else
            other_sites = sitelist[1:i-1]
        end
        total_schedule = vcat(other_sites, total_schedule)
        total_schedule = generateScheduleByDFS!(site, total_schedule)
        total_schedule = total_schedule[length(sitelist):end] # Strip other sitelist prepend
    end

    schedule = convert(Schedule, total_schedule, sumProductRule!)
    for entry in schedule
        if entry.node.interfaces[entry.outbound_interface_id] in sitelist
            entry.rule = expectationRule!
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
        for iteration_count = 1:algorithm.n_iterations
            execute(algorithm.schedule)
            # Check stopping criteria
            if algorithm.callback()
                break
            end
        end
    end

    algo = ExpectationPropagation(exec, schedule, sitelist, n_iterations, callback)
    inferDistributionTypes!(algo, recognition_distributions)

    return algo
end

############################################
# Type inference and preparation
############################################

function inferDistributionTypes!(algo::ExpectationPropagation, recognition_distributions::Dict{Interface,DataType})
    # Infer the payload types for all messages in algo.schedule
    # Fill schedule_entry.inbound_types and schedule_entry.outbound_type
    schedule_entries = Dict{Interface, ScheduleEntry}() # Lookup table from interface to schedule entry

    for entry in algo.schedule
        collectInboundTypes!(entry, schedule_entries, recognition_distributions, algo) # Fill entry.inbound_types
        inferOutboundType!(entry) # Infer the outbound message type, fill entry.outbound_type

        outbound_interface = entry.node.interfaces[entry.outbound_interface_id]
        schedule_entries[outbound_interface] = entry # Add entry to lookup table
    end

    return algo
end

function collectInboundTypes!(  entry::ScheduleEntry,
                                schedule_entries::Dict{Interface, ScheduleEntry},
                                recognition_distributions::Dict{Interface,DataType},
                                algo::ExpectationPropagation)
    # Look up the types of the inbound messages for entry.
    # Fill entry.inbound_types
    entry.inbound_types = []

    for (id, interface) in enumerate(entry.node.interfaces)
        if (id == entry.outbound_interface_id) && (entry.rule == sumProductRule!)
            # Incoming msg on outbound interface is always Void for sumProductRule! rule
            push!(entry.inbound_types, Void)
        elseif haskey(recognition_distributions, interface.partner)
            # Incoming msg from a site, so the type is given by the recognition distribution
            push!(entry.inbound_types, Message{recognition_distributions[interface.partner]})
        else
            # Incoming msg from earlier schedule entry
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

function compile!(entry::ScheduleEntry, ::Type{Val{symbol(expectationRule!)}}, ::InferenceAlgorithm)
    # Generate entry.execute for schedule entry with expectationRule! calculation rule

    inbound_messages = [interface.partner.message for interface in entry.node.interfaces]

    return buildExecute!(entry, inbound_messages)
end
