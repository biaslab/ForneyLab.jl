#####################
# Unit tests
#####################

facts("TerminalNode unit tests") do
    context("PriorNode should be a type alias to TerminalNode") do
        @fact PriorNode --> TerminalNode
    end

    context("TerminalNode() should initialize a TerminalNode with 1 interface") do
        FactorGraph()
        TerminalNode(id=:node)
        @fact length(n(:node).interfaces) --> 1
        @fact n(:node).i[:out] --> n(:node).interfaces[1]
        @fact ForneyLab.firstFreeInterface(n(:node)) --> n(:node).i[:out]
        Edge(n(:node), TerminalNode())
        @fact_throws ForneyLab.firstFreeInterface(n(:node))
    end

    context("TerminalNode should throw an error when its value is set to a message") do
        FactorGraph()
        @fact TerminalNode(DeltaDistribution(1.0)).value --> DeltaDistribution(1.0)
        @fact_throws TerminalNode(Message(DeltaDistribution(1.0)))
    end

    context("TerminalNode should propagate distributions") do
        validateOutboundMessage(TerminalNode(GaussianDistribution(m=4.0, V=5.0)),
                                1,
                                [nothing],
                                GaussianDistribution(m=4.0, V=5.0))

        validateOutboundMessage(TerminalNode(DeltaDistribution(4.0)),
                                1,
                                [nothing],
                                DeltaDistribution(4.0))
    end
end
