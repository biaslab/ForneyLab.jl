export
variationalExpectationPropagationAlgorithm

"""
Create a variational EP algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(algo::Algorithm; free_energy=false)
    (length(algo.recognition_factors) > 0) || error("No recognition factors defined on algorithm. Pass a factorization or factorized Algorithm object to create a variational algorithm.")
    for (id, rf) in algo.recognition_factors
        # Set the target regions (variables and clusters) of the recognition factor
        setTargets!(rf, algo, free_energy=free_energy, external_targets=true)

        # Infer schedule and marginal computations
        schedule = variationalExpectationPropagationSchedule(rf)
        rf.schedule = condense(flatten(schedule)) # Inline all internal message passing and remove clamp node entries
        rf.marginal_table = marginalTable(rf)
    end

    # Populate fields for algorithm compilation
    assembleAlgorithm!(algo)
    free_energy && assembleFreeEnergy!(algo)

    return algo
end

function variationalExpectationPropagationAlgorithm(args::Vararg{Union{T, Set{T}, Vector{T}} where T<:Variable}; ids=Symbol[], id=Symbol(""), free_energy=false)
    rfz = Algorithm(args...; ids=ids, id=id)
    algo = variationalExpectationPropagationAlgorithm(rfz, free_energy=free_energy)

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
    schedule = summaryPropagationSchedule(sort(collect(recognition_factor.target_variables), rev=true); target_sites=[breaker_sites; ep_sites], limit_set=internal_edges)

    nodes_connected_to_external_edges = nodesConnectedToExternalEdges(recognition_factor)
    for entry in schedule
        if entry.interface in ep_sites
            entry.message_update_rule = ExpectationPropagationRule{typeof(entry.interface.node)}
        elseif (entry.interface.node in nodes_connected_to_external_edges) && !isa(entry.interface.node, DeltaFactor)
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
