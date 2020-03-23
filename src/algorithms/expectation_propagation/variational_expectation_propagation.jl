export
variationalExpectationPropagationAlgorithm

"""
Create a variational EP algorithm to infer marginals over a posterior distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(pfz::PosteriorFactorization=currentPosteriorFactorization(); 
                                                    id=Symbol(""), 
                                                    free_energy=false)

    (length(pfz.posterior_factors) > 0) || error("No factors defined on posterior factorization.")

    # Set the target regions (variables and clusters) for each posterior factor
    for (_, pf) in pfz.posterior_factors
        setTargets!(pf, pfz, free_energy=free_energy, external_targets=true)
    end

    # Infer schedule and marginal computations for each recogition factor
    for (_, pf) in pfz.posterior_factors
        schedule = variationalExpectationPropagationSchedule(pf)
        pf.schedule = condense(flatten(schedule)) # Inline all internal message passing and remove clamp node entries
        pf.marginal_table = marginalTable(pf)
    end

    # Populate fields for algorithm compilation
    algo = InferenceAlgorithm(pfz, id=id)
    assembleInferenceAlgorithm!(algo)
    free_energy && assembleFreeEnergy!(algo)

    return algo
end

function variationalExpectationPropagationAlgorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[], id=Symbol(""), free_energy=false)
    pfz = PosteriorFactorization(args...; ids=ids, id=id)
    algo = variationalExpectationPropagationAlgorithm(pfz, id=id, free_energy=free_energy)

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

    # Schedule messages towards posterior distributions and target sites, limited to the internal edges
    schedule = summaryPropagationSchedule(sort(collect(posterior_factor.target_variables), rev=true); target_sites=[breaker_sites; ep_sites], limit_set=internal_edges)

    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(posterior_factor)
    for entry in schedule
        if entry.interface in ep_sites
            entry.message_update_rule = ExpectationPropagationRule{typeof(entry.interface.node)}
        elseif (entry.interface.node in nodes_connected_to_external_edges) && !isa(entry.interface.node, DeltaFactor)
            local_posterior_factors = localPosteriorFactors(entry.interface.node)
            if allunique(local_posterior_factors) # Local posterior factorization is naive
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
