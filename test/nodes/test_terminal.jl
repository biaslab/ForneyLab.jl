#####################
# Unit tests
#####################

facts("TerminalNode unit tests") do
    context("PriorNode should be a type alias to TerminalNode") do
        @fact PriorNode => TerminalNode
    end

    context("TerminalNode() should initialize a TerminalNode with 1 interface") do
        FactorGraph()
        node = TerminalNode()
        @fact typeof(node) => TerminalNode
        @fact length(node.interfaces) => 1
        @fact node.i[:out] => node.interfaces[1]
        @fact ForneyLab.firstFreeInterface(node) => node.i[:out]
        Edge(node, TerminalNode())
        @fact_throws ForneyLab.firstFreeInterface(node)
    end

    context("TerminalNode should throw an error when its value is set to a message") do
        @fact TerminalNode(DeltaDistribution(1.0)).value => DeltaDistribution(1.0)
        @fact_throws TerminalNode(Message(DeltaDistribution(1.0)))
    end

    context("TerminalNode should propagate a GaussianDistribution") do
        node = TerminalNode(GaussianDistribution(m=2.0, V=4.0))
        @fact node.interfaces[1].message => nothing
        (rule, msg) = ForneyLab.sumProduct!(node, 1, nothing)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) => Message{GaussianDistribution}
        @fact node.interfaces[1].message.payload.m => [2.0]
        @fact node.interfaces[1].message.payload.V => reshape([4.0], 1, 1)
    end

    context("TerminalNode should propagate a DeltaDistribution") do
        node = TerminalNode(DeltaDistribution([1.0, 2.0]))
        @fact node.interfaces[1].message => nothing
        (rule, msg) = ForneyLab.sumProduct!(node, 1, nothing)
        @fact node.interfaces[1].message => msg
        @fact typeof(node.interfaces[1].message) <: Message => true
        @fact typeof(node.interfaces[1].message.payload) <: DeltaDistribution => true
        @fact node.interfaces[1].message.payload.m => [1.0, 2.0]
    end
end