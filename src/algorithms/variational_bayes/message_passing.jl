function execute(subgraph::Subgraph, q_distributions::Dict{Tuple{Node, Subgraph}, QDistribution}, factorization::QFactorization)
    # Execute internal schedule
    execute(subgraph.internal_schedule)

    # Update q-distributions at external edges
    if ForneyLab.verbose
        println(" ")
        println("Marginals (node), result")
        println("--------------------------------------------")
    end
    for node in subgraph.external_schedule
        d = calculateQDistribution!(q_distributions, node, subgraph, factorization)
        if ForneyLab.verbose
            println("$(node.id) : $(format(d))")
        end
    end
end

function execute(factorization::QFactorization, q_distributions::Dict{Tuple{Node, Subgraph}, QDistribution})
    for subgraph in factorization.factors
        execute(subgraph, q_distributions, factorization)
    end
end
