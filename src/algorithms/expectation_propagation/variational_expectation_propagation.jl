export
variationalExpectationPropagationAlgorithm,
variationalExpectationPropagationSchedule

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(pfz::PosteriorFactorization=currentPosteriorFactorization(), 
    id=Symbol(""))
    
    for (id, pf) in pfz
        schedule = variationalExpectationPropagationSchedule(pf)
        pf.schedule = condense(flatten(schedule)) # Inline all internal message passing and remove clamp node entries
        pf.marginal_table = marginalTable(pf)
    end

    algo = InferenceAlgorithm(pfz,id)
    assembleInferenceAlgorithm!(algo)

    return algo
end

function variationalExpectationPropagationAlgorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[], id=Symbol(""))
    pfz = PosteriorFactorization(args...; ids=ids)
    algo = variationalExpectationPropagationAlgorithm(pfz, id)

    return algo
end

"""
`variationalExpectationPropagationSchedule()` generates an expectation propagation message passing schedule
that is limited to the `posterior_factor`. Updates on EP sites are computed with an `ExpectationPropagationRule`.
"""
function variationalExpectationPropagationSchedule(posterior_factor::PosteriorFactor)
    internal_edges = posterior_factor.internal_edges
    ep_sites = collectEPSites(nodes(internal_edges))
    breaker_sites = Interface[site.partner for site in ep_sites]
    breaker_types = breakerTypes(breaker_sites)

    # Schedule messages towards recognition distributions and target sites, limited to the internal edges
    schedule = summaryPropagationSchedule(sort(collect(posterior_factor.variables), rev=true); target_sites=[breaker_sites; ep_sites], limit_set=internal_edges)

    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(posterior_factor)
    for entry in schedule
        if entry.interface in ep_sites
            entry.message_update_rule = ExpectationPropagationRule{typeof(entry.interface.node)}
        elseif entry.interface.node in nodes_connected_to_external_edges
            local_posterior_factors = localPosteriorFactors(entry.interface.node)
            if allunique(local_posterior_factors) # Local recognition factorization is naive
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
