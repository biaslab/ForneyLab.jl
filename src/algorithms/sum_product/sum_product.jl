export SumProduct

include("scheduler.jl")

type SumProduct <: InferenceAlgorithm
    execute::Function
    schedule::Schedule
end

#--------------------------------
# SumProduct constructors
#--------------------------------

function SumProduct(graph::FactorGraph=currentGraph())
    # Generates a SumProduct algorithm that propagates messages to all wraps and write buffers.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(graph)
    exec(algorithm) = execute(algorithm.schedule)
    return SumProduct(exec, schedule)
end

function SumProduct(outbound_interface::Interface)
    # Generates a SumProduct algorithm to calculate the outbound message on outbound_interface.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(outbound_interface)
    exec(algorithm) = execute(algorithm.schedule)
    return SumProduct(exec, schedule)
end

function SumProduct(partial_list::Vector{Interface})
    # Generates a SumProduct algorithm that at least propagates to all interfaces in the argument vector.
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule(partial_list)
    exec(algorithm) = execute(algorithm.schedule)
    return SumProduct(exec, schedule)
end

function SumProduct(edge::Edge)
    # Generates a SumProduct algorithm to calculate the marginal on edge
    # Only works in acyclic graphs.
    schedule = generateSumProductSchedule([edge.head, edge.tail])
    function exec(algorithm)
        execute(algorithm.schedule)
        calculateMarginal!(edge)
    end
    return SumProduct(exec, schedule)
end


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface, ::Type{Val{symbol("ForneyLab.sumProduct!")}})
    # Sum-product specific method to collect all required inbound messages in an array.
    # This array is used to call the node update function (sumProduct!).
    # outbound_interface: the interface on which the outbound message will be updated.
    # If include_inbound_on_outbound_interface is true, the inbound message on the outbound interface will also be required and included in the array.
    # Returns: (outbound interface id, array of inbound messages).

    outbound_interface_index = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            outbound_interface_index = j
            push!(inbounds, nothing)
            continue
        end

        try
            push!(inbounds, interface.partner.message)
        catch
            error("Cannot collect inbound message on $(interface). Make sure there is an inbound message present at this interface.")
        end
    end

    return (outbound_interface_index, inbounds)
end
