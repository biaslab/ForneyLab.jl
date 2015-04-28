function execute(subgraph::Subgraph, q_distributions::Dict{(Node, Subgraph), QDistribution}, factorization::QFactorization)
    ForneyLab.printVerbose("Subgraph $(subgraph):")

    # Execute internal schedule
    SumProduct.execute(subgraph.internal_schedule)

    # Update q-distributions at external edges
    for node in subgraph.external_schedule
        calculateQDistribution!(q_distributions, node, subgraph, factorization)
    end
end

function execute(factorization::QFactorization, q_distributions::Dict{(Node, Subgraph), QDistribution})
    for subgraph in factorization.factors
        execute(subgraph, q_distributions, factorization)
    end
end
