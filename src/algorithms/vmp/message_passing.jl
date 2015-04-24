function execute(subgraph::Subgraph, q_distributions::Dict{(Node, Subgraph), ProbabilityDistribution})
    printVerbose("Subgraph $(findfirst(scheme.factorization, subgraph)):")

    # Execute internal schedule
    execute(subgraph.internal_schedule)

    # Update q-distributions at external edges
    for node in subgraph.external_schedule
        calculateQDistribution!(q_distributions, node, subgraph)
    end
end

function execute(factorization::Factorization, q_distributions::Dict{(Node, Subgraph), ProbabilityDistribution})
    for subgraph in factorization.factors
        execute(subgraph, q_distributions)
    end
end
