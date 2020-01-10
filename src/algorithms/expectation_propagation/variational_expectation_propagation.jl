export
variationalExpectationPropagationAlgorithm,
variationalExpectationPropagationSchedule

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(algo::Algorithm=currentAlgorithm())
    for (id, rf) in algo.recognition_factors
        schedule = variationalExpectationPropagationSchedule(rf)
        rf.schedule = condense(flatten(schedule)) # Inline all internal message passing and remove clamp node entries
        rf.marginal_table = marginalTable(rf)
    end

    assembleAlgorithm!(algo)

    return algo
end

function variationalExpectationPropagationAlgorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[])
    rfz = Algorithm(args..., ids=ids)
    algo = variationalExpectationPropagationAlgorithm(rfz)

    return algo
end

"""
`variationalExpectationPropagationSchedule()` generates an expectation propagation message passing schedule
that is limited to the `recognition_factor`. Updates on EP sites are computed with an `ExpectationPropagationRule`.
"""
function variationalExpectationPropagationSchedule(recognition_factor::RecognitionFactor)
    internal_edges = recognition_factor.internal_edges
    ep_sites = collectEPSites(nodes(internal_edges))
    breaker_sites = Interface[site.partner for site in ep_sites]
    breaker_types = breakerTypes(breaker_sites)

    # Schedule messages towards recognition distributions and target sites, limited to the internal edges
    schedule = summaryPropagationSchedule(sort(collect(recognition_factor.variables), rev=true); target_sites=[breaker_sites; ep_sites], limit_set=internal_edges)

    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(recognition_factor)
    for entry in schedule
        if entry.interface in ep_sites
            entry.message_update_rule = ExpectationPropagationRule{typeof(entry.interface.node)}
        elseif entry.interface.node in nodes_connected_to_external_edges
            local_recognition_factors = localRecognitionFactors(entry.interface.node)
            if allunique(local_recognition_factors) # Local recognition factorization is naive
                entry.message_update_rule = NaiveVariationalRule{typeof(entry.interface.node)}
            else
                entry.message_update_rule = StructuredVariationalRule{typeof(entry.interface.node)}
            end
        else
            entry.message_update_rule = SumProductRule{typeof(entry.interface.node)}
        end
    end

    inferUpdateRules!(schedule, inferred_outbound_types=breaker_types)

    return schedule
end
