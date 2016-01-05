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
            callback::Function = () -> false)
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

    # Build schedule for step 1
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

    function exec(algorithm)
        # Init all sites with vague messages
        for site in algorithm.sites
            site.message = Message(vague!(site.edge.distribution_type))
        end
        for iteration_count = 1:algorithm.num_iterations
            execute(algorithm.schedule)
            # Check stopping criteria
            if algorithm.callback()
                break
            end
        end
    end

    algo = ExpectationPropagation(exec, schedule, sites, num_iterations, callback)
    compile!(algo.schedule, algo)

    return algo
end

function compile!(schedule_entry::ScheduleEntry, ::Type{Val{symbol("ForneyLab.ep!")}}, ::InferenceAlgorithm)
    # Compile ScheduleEntry objects for SumProduct algorithm
    # Generates schedule_entry.execute function

    # Collect references to all required inbound messages for executing message computation rule
    outbound_interface = schedule_entry.interface
    node = outbound_interface.node
    outbound_interface_index = 0
    inbounds = Message[]
    for j = 1:length(node.interfaces)
        interface = node.interfaces[j]
        if is(interface, outbound_interface)
            outbound_interface_index = j
        end
        push!(inbounds, interface.partner.message)
    end

    # Assign the "compiled" update rule as an anomynous function to the schedule entry execute field
    schedule_entry.execute = ( () -> ep!(node, outbound_interface_index, inbounds...) )

    return schedule_entry
end
