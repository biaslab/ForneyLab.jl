module ExpectationPropagation

using ..ForneyLab

#--------------------------------
# Algorithm constructors
#--------------------------------

# function Algorithm(graph::FactorGraph=currentGraph())
#
# end

# function Algorithm(outbound_interface::Interface)
#
# end
#
# function Algorithm(edge::Edge)
#
# end

function Algorithm(
            sites::Vector{Interface};
            num_iterations::Int64 = 100000,
            callback::Function = () -> false)
    # Build an EP message passing algorithm for the specified sites
    # num_iterations specifies the maximum number of iterations
    # After each iteration, callback is called
    # If this function returns true, the algorithm stops

    # Algorithm overview:
    # 1. Calculate all required messages that do not depend on the sites
    # 2. Init all sites with vague messages
    # 3. Calculate all cavity distributions
    # 4. Approximate all site distributions
    # 5. Check stopping criteria, goto 3

    (length(sites) > 0) || error("Specify at least one site")

    # Build schedule for step 1
    pre_sites_interfaces = Vector{Interface}()
    for site in sites
        pre_site_interfaces = SumProduct.generateScheduleByDFS!(site)
        for idx=1:length(pre_site_interfaces)
            if pre_site_interfaces[idx] in sites
                pre_sites_interfaces = vcat(pre_sites_interfaces, pre_site_interfaces[1:idx-1])
                break
            end
        end
    end

    # Build schedule for step 3
    cavity_interfaces = [site.partner for site in sites]
    cavities_dependencies = vcat(pre_sites_interfaces, sites) # init with schedule that has already been executed once we reach step 3
    pre_cavities_schedule_length = length(cavities_dependencies)
    for cavity_interface in [site.partner for site in sites]
        cavities_dependencies = SumProduct.generateScheduleByDFS!(cavity_interface, cavities_dependencies)
    end
    if pre_cavities_schedule_length < length(cavities_dependencies)
        # Get rid of the init prefix from step 1 and 2
        cavities_dependencies = cavities_dependencies[pre_cavities_schedule_length+1:end]
    else
        error("Cannot construct a schedule for the cavity distributions")
    end

    # Build schedule for sites
    ep_dependencies = vcat(pre_sites_interfaces, cavities_dependencies)
    pre_ep_schedule_length = length(ep_dependencies)
    for site in sites
        ep_dependencies = SumProduct.generateScheduleByDFS!(site, ep_dependencies)
    end
    sites_schedule = convert(Schedule, ep_dependencies[pre_ep_schedule_length+1:end], ep!)
    for ep_schedule_entry in sites_schedule
        if !(ep_schedule_entry.interface in sites)
            ep_schedule_entry.message_calculation_rule = sumProduct!
        end
    end

    fields = Dict{Symbol,Any}(
                :sites => sites,
                :pre_sites_schedule => convert(Schedule, pre_sites_interfaces, sumProduct!),
                :cavities_schedule => convert(Schedule, cavities_dependencies, sumProduct!),
                :sites_schedule => sites_schedule,
                :num_iterations => num_iterations,
                :callback => callback
                )

    function exec(fields)
        # 1. Calculate all required messages that do not depend on the sites
        execute(fields[:pre_sites_schedule])
        # 2. Init all sites with vague messages
        for site in fields[:sites]
            site.message = Message(vague(site.edge.distribution_type))
        end

        for iteration_count = 1:fields[:num_iterations]
            # 3. Calculate all cavity distributions
            execute(fields[:cavities_schedule])
            # 4. Approximate all site distributions
            execute(fields[:sites_schedule])
            # 5. Check stopping criteria
            if fields[:callback]()
                break
            end
        end

    end

    return ForneyLab.Algorithm(exec, fields)
end


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface)
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

end # module
