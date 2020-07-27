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
        schedule = messagePassingSchedule(pf)
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