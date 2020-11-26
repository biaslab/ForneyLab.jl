export InferenceAlgorithm, currentInferenceAlgorithm, messagePassingAlgorithm

"""
An `InferenceAlgorithm` specifies the computations for the quantities of interest.
"""
mutable struct InferenceAlgorithm
    id::Symbol
    posterior_factorization::PosteriorFactorization

    # Bookkeeping for faster lookup during assembly
    interface_to_schedule_entry::Dict{Interface, ScheduleEntry}
    target_to_marginal_entry::Dict{Union{Variable, Cluster}, MarginalEntry}

    # Fields for free energy algorithm assembly
    average_energies::Vector{Dict{Symbol, Any}}
    entropies::Vector{Dict{Symbol, Any}}
end

"""
Return currently active `InferenceAlgorithm`.
Create one if there is none.
"""
function currentInferenceAlgorithm()
    try
        return current_inference_algorithm
    catch
        return InferenceAlgorithm()
    end
end

function setCurrentInferenceAlgorithm(algo::InferenceAlgorithm)
     global current_inference_algorithm = algo
end

function InferenceAlgorithm(
    pfz=currentPosteriorFactorization();
    id=Symbol(""))

    setCurrentInferenceAlgorithm(
        InferenceAlgorithm(
            id,
            pfz,
            Dict{Interface, ScheduleEntry}(),
            Dict{Union{Variable, Cluster}, MarginalEntry}(),
            Dict{Symbol, Any}[],
            Dict{Symbol, Any}[]))
end

"""
Create a message passing algorithm to infer marginals over a posterior distribution
"""
function messagePassingAlgorithm(target_variables::Vector{Variable}=Variable[], # Quantities of interest
                                 pfz::PosteriorFactorization=currentPosteriorFactorization(); 
                                 ep_sites=Tuple[],
                                 id=Symbol(""), 
                                 free_energy=false)

    if isempty(pfz.posterior_factors) # If no factorization is defined
        PosteriorFactor(pfz.graph, pfz=pfz, id=Symbol("")) # Contain the entire graph in a single posterior factor
    end

    # Set the EP sites for each corresponding posterior factor
    for ep_site in ep_sites
        site_node = pfz.graph.nodes[ep_site[1]] # Find node by id
        site_iface = site_node.i[ep_site[2]] # Find interface by id
        site_edge = site_iface.edge
        for (_, pf) in pfz.posterior_factors # Find corresponding posterior factor
            if site_edge in pf.internal_edges
                push!(pf.ep_sites, site_iface)
                break # Skip to next site
            end
        end
    end

    # Set the targets for each posterior factor
    for (_, pf) in pfz.posterior_factors
        setTargets!(pf, pfz, target_variables=Set(target_variables), free_energy=free_energy, external_targets=true)
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

messagePassingAlgorithm(target_variable::Variable,
                        pfz::PosteriorFactorization=currentPosteriorFactorization(); 
                        ep_sites=Tuple[],
                        id=Symbol(""), 
                        free_energy=false) = messagePassingAlgorithm([target_variable], pfz; ep_sites=ep_sites, id=id, free_energy=free_energy)

function interfaceToScheduleEntry(algo::InferenceAlgorithm)
    mapping = Dict{Interface, ScheduleEntry}()
    for (id, pf) in algo.posterior_factorization
        pf_mapping = interfaceToScheduleEntry(pf.schedule)
        merge!(mapping, pf_mapping)
    end

    return mapping
end

function targetToMarginalEntry(algo::InferenceAlgorithm)
    mapping = Dict{Union{Cluster, Variable}, MarginalEntry}()
    for (id, pf) in algo.posterior_factorization
        pf_mapping = targetToMarginalEntry(pf.marginal_table)
        merge!(mapping, pf_mapping)
    end

    return mapping    
end