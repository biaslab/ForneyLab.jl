function execute(subgraph::Subgraph, recognition_distributions::Dict{Tuple{Node, Subgraph}, RecognitionDistribution}, factorization::RecognitionFactorization)
    # Execute internal schedule
    execute(subgraph.internal_schedule)

    # Update recognition distributions at external edges
    if ForneyLab.verbose
        println(" ")
        println("Marginals (node), result")
        println("--------------------------------------------")
    end
    for node in subgraph.external_schedule
        d = calculateRecognitionDistribution!(recognition_distributions, node, subgraph, factorization)
        if ForneyLab.verbose
            println("$(node.id) : $(format(d))")
        end
    end
end

function execute(factorization::RecognitionFactorization, recognition_distributions::Dict{Tuple{Node, Subgraph}, RecognitionDistribution})
    for subgraph in factorization.factors
        execute(subgraph, recognition_distributions, factorization)
    end
end
