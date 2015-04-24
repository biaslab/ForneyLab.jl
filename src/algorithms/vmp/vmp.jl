module VMP

using ..ForneyLab

include("generate_schedule.jl")
include("subgraph.jl")
include("factorization.jl")
include("generate_schedule.jl")
#include("calculate_q_distribution.jl")
#include("message_passing.jl")

#--------------------------------
# Algorithm specific constructors
#--------------------------------

function Algorithm(graph::FactorGraph=currentGraph())
    # Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on graph
    # Uses a mean field factorization and autoscheduler
    factorization = factorize(graph) # Mean field factorization
    generateSchedule!(factorization, graph) # Generate and store internal and external schedules on factorization subgraphs
    q_distributions = vagueQDistributions(factorization) # Initialize vague q distributions

    # Construct the algorithm execute function
    exec(fields) = execute(factorization, q_distributions, graph) # For all subgraphs, execute internal and external schedules
    return ForneyLab.Algorithm(exec, {:factorization => factorization, :q_distributions => q_distributions})
end


#---------------------------------------------------
# Construct algorithm specific update-call signature
#---------------------------------------------------

function collectInbounds(outbound_interface::Interface, algorithm::Algorithm=current_algorithm)
    # VMP specific method to collect all required inbound messages and marginals in an array.
    # This array is used to call the node update function (vmp!)
    # outbound_interface: the interface on which the outbound message will be updated
    # Returns: (outbound interface id, array of inbound messages and marginals)
    
    outbound_interface_id = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = outbound_interface.node.interfaces[j]
        if is(interface, outbound_interface)
            # We don't need the inbound message on the outbound interface
            outbound_interface_id = j
            push!(inbounds, nothing) # This interface is outbound, push "nothing"
        else
            if !haskey(algorithm.fields[:factorization].edge_to_subgraph, interface.edge) || !haskey(algorithm.fields[:factorization].edge_to_subgraph, outbound_interface.edge)
                # Inbound and/or outbound edge is not explicitly listed in the algorithm.fields.
                # This is possible if one of those edges is internal to a composite node.
                # We will default to sum-product message passing, and consume the message on the inbound interface.
                # Composite nodes with explicit message passing will throw an error when one of their external interfaces belongs to a different subgraph, so it is safe to assume sum-product.
                try return push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            end

            # Should we require the inbound message or marginal?
            if is(algorithm.fields[:factorization].edge_to_subgraph[interface.edge], algorithm.fields[:factorization].edge_to_subgraph[outbound_interface.edge])
                # Both edges in same subgraph, require message
                try push!(inbounds, interface.partner.message) catch error("$(interface) is not connected to an edge.") end
            else
                # A subgraph border is crossed, require marginal
                # The factor is the set of internal edges that are in the same subgraph
                sg = algorithm.fields[:factorization].edge_to_subgraph[interface.edge]
                try push!(inbounds, algorithm.fields[:q_distributions][(node, sg)]) catch error("Missing approximate marginal for $(interface)") end
            end
        end
    end

    return (outbound_interface_id, inbounds)
end

end