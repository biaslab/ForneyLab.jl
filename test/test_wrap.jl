facts("Wrap integration tests") do
    # Wrap()
    context("Wrap() should register a wrap for a pair of TerminalNodes") do
        g = initializeBufferGraph()
        wrap = Wrap(n(:node_t2), n(:node_t1))
        @fact length(g.wraps) --> 1
        @fact wrap in values(g.wraps) --> true
        @fact wrap.id --> :node_t2_node_t1
        @fact wrap.tail --> n(:node_t2)
        @fact wrap.head --> n(:node_t1)
        @fact g.current_section --> 1
        @fact length(wrap.head_buffer) --> 1
        @fact length(wrap.tail_buffer) --> 1
    end

    context("Wrap() should validate and accept an id") do
        g = initializeWrapGraph()
        wrap = Wrap(n(:t2), n(:t1), id=:my_wrap)
        @fact wrap.id --> :my_wrap
        @fact_throws Wrap(n(:t4), n(:t3), id=:my_wrap) # Same name
        @fact_throws Wrap(n(:t4), n(:t1)) # Multiple sinks
    end

    # wrap
    context("wrap() should find a wrap by id") do
        g = initializeBufferGraph()
        wr = Wrap(n(:node_t2), n(:node_t1), id=:my_wrap)
        @fact wrap(:my_wrap) --> wr
    end

    # wraps
    context("wraps() should return wraps") do
        g = initializeWrapGraph()
        wrap1 = Wrap(n(:t2), n(:t1))
        wrap2 = Wrap(n(:t2), n(:t3))
        @fact wraps(g) --> Set{Wrap}([wrap1, wrap2])
        @fact wraps(n(:t1)) --> Set{Wrap}([wrap1])
        @fact wraps(n(:t2)) --> Set{Wrap}([wrap1, wrap2])
    end

    context("Wrap() with specified block_size should initialize corresponding field on the graph and tail/head buffers") do
        g = initializeWrapGraph()
        wrap1 = Wrap(n(:t2), n(:t1), block_size=10)
        @fact isdefined(g, :block_size) --> true
        @fact g.block_size --> 10
        @fact length(wrap1.head_buffer) --> g.block_size
        @fact length(wrap1.tail_buffer) --> g.block_size
    end
    
end
