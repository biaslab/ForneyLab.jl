module SumProduct

using ..ForneyLab

include("generate_schedule.jl")


#--------------------------------
# Algorithm specific constructors
#--------------------------------

function Algorithm(graph::FactorGraph=currentGraph())
    # Generates a sumproduct algorithm
    # Uses autoscheduler and only works in acyclic graphs
    schedule = SumProduct.generateSchedule(graph)

    # Construct the execute function and its arguments
    exec(fields) = execute(fields[:schedule])
    return ForneyLab.Algorithm(exec, {:schedule => schedule})
end

function Algorithm(outbound_interface::Interface)
    # Generates a sumproduct algorithm to calculate the outbound message on outbound_interface
    # Uses autoscheduler and only works in acyclic graphs
    schedule = SumProduct.generateSchedule(outbound_interface)

    # Construct the execute function and its arguments
    exec(fields) = execute(fields[:schedule])
    return ForneyLab.Algorithm(exec, {:schedule => schedule})
end

function Algorithm(partial_list::Vector{Interface})
    # Generates a sumproduct algorithm that at least propagates to all interfaces in the argument vector.
    # Uses autoscheduler and only works in acyclic graphs
    schedule = SumProduct.generateSchedule(partial_list)

    # Construct the execute function and its arguments
    exec(fields) = execute(fields[:schedule])
    return ForneyLab.Algorithm(exec, {:schedule => schedule})
end

function Algorithm(edge::Edge)
    # Generates a sumproduct algorithm to calculate the marginal on edge
    # Uses autoscheduler and only works in acyclic graphs
    schedule = SumProduct.generateSchedule([edge.head, edge.tail])

    # Construct the execute function and its arguments
    function exec(fields)
        execute(fields[:schedule])
        calculateMarginal!(fields[:edge])
    end
    return ForneyLab.Algorithm(exec, {:schedule => schedule, :edge => edge})
end


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface)
    # Sum-product specific method to collect all required inbound messages in an array.
    # This array is used to call the node update function (sumProduct!)
    # outbound_interface: the interface on which the outbound message will be updated
    # Returns: (outbound interface id, array of inbound messages)
    
    outbound_interface_index = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            # We don't need the inbound message on the outbound interface
            outbound_interface_index = j
            push!(inbounds, nothing) # This interface is outbound, push "nothing"
        else
            try 
                push!(inbounds, interface.partner.message) 
            catch 
                error("Cannot collect inbound message on $(interface). Make sure there is an inbound message present at this interface.") 
            end
        end
    end

    return (outbound_interface_index, inbounds)
end

end