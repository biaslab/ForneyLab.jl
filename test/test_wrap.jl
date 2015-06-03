facts("Wrap integration tests") do
    # Wrap()
    context("Wrap() should register a timewrap for a pair of TerminalNodes") do
        g = initializeBufferGraph()
        wrap = Wrap(n(:node_t1), n(:node_t2))
        @fact length(g.wraps) => 1
        @fact wrap in g.wraps => true
    end

    # clearWraps
    context("clearWraps should deregister all time wraps") do
        g = initializeBufferGraph()
        Wrap(n(:node_t1), n(:node_t2))
        clearWraps(g)
        @fact length(g.wraps) => 0
    end
end