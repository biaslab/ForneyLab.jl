#####################
# Unit tests
#####################

facts("BernoulliNode unit tests") do
    context("BernoulliNode should initialize a BernoulliNode with 2 interfaces") do
        FactorGraph()
        BernoulliNode(id=:node)
        @fact length(n(:node).interfaces) --> 2
        @fact n(:node).i[:in] --> n(:node).interfaces[1]
        @fact n(:node).i[:out] --> n(:node).interfaces[2]
    end

    context("U() should evaluate the average energy") do
        @fact ForneyLab.U(BernoulliNode, Beta(a=1.0, b=2.0), Bernoulli(0.2)) --> roughly(0.7)
    end    

    context("BernoulliNode should pass sum-product messages") do
        # Forward message
        validateOutboundMessage(BernoulliNode(),
                                2,
                                [Message(Beta(a=1.0, b=2.0)), nothing],
                                Bernoulli(1/3))
        validateOutboundMessage(BernoulliNode(),
                                2,
                                [Message(Delta(0.2)), nothing],
                                Bernoulli(0.2))
        # Backward message
        validateOutboundMessage(BernoulliNode(),
                                1,
                                [nothing, Message(Delta(true))],
                                Beta(a=2.0, b=1.0))
    end

    context("BernoulliNode should pass variational messages") do
        # Forward message
        validateOutboundMessage(BernoulliNode(),
                                2,
                                [Beta(a=1.0, b=2.0), nothing],
                                Bernoulli(0.268941421369995),
                                ForneyLab.variationalRule!)
        # Backward message
        validateOutboundMessage(BernoulliNode(),
                                1,
                                [nothing, Bernoulli(0.8)],
                                Beta(a=1.8, b=1.2),
                                ForneyLab.variationalRule!)
    end
end


#####################
# Integration tests
#####################

facts("BernoulliNode integration tests") do
    context("Sum-product message passing should estimate Beta distribution parameters") do
        data = [true, false, false, true, true]

        FactorGraph()

        TerminalNode(Beta(a=0.0, b=0.0), id=:prior)
        TerminalNode(Beta(a=0.0, b=0.0), id=:term)
        EqualityNode(id=:eq)
        BernoulliNode(id=:bern)
        TerminalNode(Delta(true), id=:obs)

        Edge(n(:prior).i[:out], n(:eq).i[1])
        Edge(n(:eq).i[2], n(:term).i[:out])
        Edge(n(:eq).i[3], n(:bern).i[:in])
        Edge(n(:bern).i[:out], n(:obs).i[:out])

        Wrap(n(:term), n(:prior))

        attachReadBuffer(n(:obs), data)
        buff = attachWriteBuffer(n(:term).i[:out].partner)

        algo = SumProduct()
        run(algo)

        @fact buff[end] --> Beta(a=3.0, b=2.0)
    end

    context("Variational message passing should estimate Beta distribution parameters") do
        data = [Bernoulli(0.9), Bernoulli(0.1), Bernoulli(0.1), Bernoulli(0.9), Bernoulli(0.1)]
        K = length(data)

        FactorGraph()

        for k = 1:K
            EqualityNode(id=:eq_*k)
            BernoulliNode(id=:bern_*k)
            TerminalNode(data[k], id=:obs_*k)

            Edge(n(:eq_*k).i[3], n(:bern_*k).i[:in], id=:x_*k)
            Edge(n(:bern_*k).i[:out], n(:obs_*k).i[:out], id=:c_*k)

            if k > 1
                Edge(n(:eq_*(k-1)).i[2], n(:eq_*k).i[1])
            end
        end

        TerminalNode(Beta(a=1.0, b=1.0), id=:prior)
        TerminalNode(vague(Beta), id=:term)

        Edge(n(:prior).i[:out], n(:eq_1).i[1])
        Edge(n(:eq_*K).i[2], n(:term).i[:out])

        buff = attachWriteBuffer(n(:term).i[:out].partner)

        RecognitionFactorization()

        factor(eg(:x_1))
        for k = 1:K
            factor(eg(:c_*k))

            initialize(eg(:c_*k), vague(Bernoulli))
            initialize(eg(:x_*k), vague(Beta))
        end

        algo = VariationalBayes(n_iterations=50)
        run(algo)

        @fact buff[end].a --> 2.7700799520862915
        @fact buff[end].b --> 4.2299200479137085
    end
end