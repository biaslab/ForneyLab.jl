#####################
# Unit tests
#####################

facts("CategoricalNode unit tests") do
    context("CategoricalNode should initialize a CategoricalNode with 2 interfaces") do
        FactorGraph()
        CategoricalNode(id=:node)
        @fact length(n(:node).interfaces) --> 2
        @fact n(:node).i[:pi] --> n(:node).interfaces[1]
        @fact n(:node).i[:z] --> n(:node).interfaces[2]
    end

    context("averageEnergy() should evaluate the average energy") do
        @fact ForneyLab.averageEnergy(CategoricalNode, Dirichlet([1.,2.,1.]), Categorical([0.2,0.2,0.6]))  --> roughly(1.6333333333333335)
    end


    context("CategoricalNode should pass variational messages") do
        # Forward message
        validateOutboundMessage(CategoricalNode(),
                                2,
                                [Dirichlet([1.,2.,2.]), nothing],
                                Categorical([0.15536240349696354,0.42231879825151825,0.42231879825151825]),
                                ForneyLab.variationalRule!)
        #Backward message
        validateOutboundMessage(CategoricalNode(),
                                1,
                                [nothing, Categorical([0.1,0.2,0.7])],
                                Dirichlet([1.1,1.2,1.7]),
                                ForneyLab.variationalRule!)
    end
end

#####################
# Integration tests
#####################
facts("CategoricalNode integration tests") do
    context("Variational message passing should estimate Dirichlet distribution parameters") do
        data = [Categorical([0.1,0.1,0.8]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1]), Categorical([0.8,0.1,0.1])]
        K = length(data)

        FactorGraph()

        for k = 1:K
            EqualityNode(id=:eq_*k)
            CategoricalNode(id=:bern_*k)
            TerminalNode(data[k], id=:obs_*k)

            Edge(n(:eq_*k).i[3], n(:bern_*k).i[:pi], id=:x_*k)
            Edge(n(:bern_*k).i[:z], n(:obs_*k).i[:out], id=:c_*k)

            if k > 1
                Edge(n(:eq_*(k-1)).i[2], n(:eq_*k).i[1])
            end
        end

        TerminalNode(Dirichlet([1.,1.,1.]), id=:prior)
        TerminalNode(vague(Dirichlet{3}), id=:term)

        Edge(n(:prior).i[:out], n(:eq_1).i[1])
        Edge(n(:eq_*K).i[2], n(:term).i[:out])

        buff = attachWriteBuffer(n(:term).i[:out].partner)

        RecognitionFactorization()

        factor(eg(:x_1))
        for k = 1:K
            factor(eg(:c_*k))

            initialize(eg(:c_*k), vague(Categorical{3}))
            initialize(eg(:x_*k), vague(Dirichlet{3}))
        end

        algo = VariationalBayes(n_iterations=50)
        run(algo)

        @fact buff[end].alpha --> [6.,1.,1.]
    end
end
