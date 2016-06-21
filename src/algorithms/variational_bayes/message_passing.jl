function execute(subgraph::Subgraph, rf::RecognitionFactorization)
    # Execute internal schedule
    execute(subgraph.internal_schedule)

    # Update recognition distributions at external edges
    if ForneyLab.verbose
        println(" ")
        println("Marginals (node), result")
        println("--------------------------------------------")
    end
    for node in subgraph.nodes_connected_to_external_edges
        d = calculateRecognitionDistribution!(rf, node, subgraph)
        if ForneyLab.verbose
            println("$(node.id) : $(format(d))")
        end
    end
end

function execute(rf::RecognitionFactorization)
    for subgraph in rf.subgraphs
        execute(subgraph, recognition_distributions, rf)
    end
end
