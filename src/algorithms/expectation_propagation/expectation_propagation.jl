export ExpectationPropagation

type ExpectationPropagation <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
    sites::Vector{Interface}
    num_iterations::Int64
    callback::Function
end

#------------------------------------
# ExpectationPropagation constructors
#------------------------------------

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
            schedule_entry.message_calculation_rule = ep!
        end
    end

    function exec(algorithm)
        # Init all sites with vague messages
        for site in algorithm.sites
            site.message = Message(vague(site.edge.distribution_type))
        end
        for iteration_count = 1:algorithm.num_iterations
            execute(algorithm.schedule)
            # Check stopping criteria
            if algorithm.callback()
                break
            end
        end
    end

    return ExpectationPropagation(exec, schedule, sites, num_iterations, callback)
end


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface, ::Type{Val{symbol("ForneyLab.ep!")}})
    # EP specific method to collect all required inbound messages in an array.
    # This array is used to call the node update function (ep!).
    # outbound_interface: the interface on which the outbound message will be updated.
    # Returns: (outbound interface id, array of inbound messages).

    outbound_interface_index = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            outbound_interface_index = j
        end

        try
            push!(inbounds, interface.partner.message)
        catch
            error("Cannot collect inbound message on $(interface). Make sure there is an inbound message present at this interface.")
        end
    end

    return (outbound_interface_index, inbounds)
end
