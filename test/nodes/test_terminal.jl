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
        @fact typeof(n(:node)) --> TerminalNode
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

    context("TerminalNode should propagate a GaussianDistribution") do
        FactorGraph()
        TerminalNode(GaussianDistribution(m=2.0, V=4.0), id=:node)
        @fact n(:node).interfaces[1].message --> nothing
        (rule, msg) = ForneyLab.sumProduct!(n(:node), 1, nothing)
        @fact n(:node).interfaces[1].message --> msg
        @fact typeof(n(:node).interfaces[1].message) --> Message{GaussianDistribution}
        @fact n(:node).interfaces[1].message.payload.m --> 2.0
        @fact n(:node).interfaces[1].message.payload.V --> 4.0
    end

    context("TerminalNode should propagate a DeltaDistribution") do
        FactorGraph()
        TerminalNode(DeltaDistribution(2.0), id=:node)
        @fact n(:node).interfaces[1].message --> nothing
        (rule, msg) = ForneyLab.sumProduct!(n(:node), 1, nothing)
        @fact n(:node).interfaces[1].message --> msg
        @fact typeof(n(:node).interfaces[1].message) <: Message --> true
        @fact typeof(n(:node).interfaces[1].message.payload) <: DeltaDistribution --> true
        @fact n(:node).interfaces[1].message.payload.m --> 2.0
    end
end
