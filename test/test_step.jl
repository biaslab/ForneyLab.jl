#####################
# Integration tests
#####################

facts("Read/write buffer integration tests") do
    # Background
    g = FactorGraph()
    node_t1 = TerminalNode()
    node_t2 = TerminalNode()
    e = Edge(node_t1, node_t2)

    # setReadBuffer
    context("setReadBuffer should register a read buffer for a TerminalNode") do
        read_buffer = zeros(10)
        setReadBuffer(node_t1, read_buffer)
        @fact g.read_buffers[node_t1] => read_buffer
    end
    context("setReadBuffer should not register a read buffer if the specified node is not in the specified graph") do
        @fact_throws setReadBuffer(TerminalNode(), zeros(10))
    end

    # setWriteBuffer
    context("setWriteBuffer should register a write buffer for an Interface") do
        write_buffer = Array(Any, 0)
        setWriteBuffer(node_t1.out, write_buffer)
        @fact g.write_buffers[node_t1.out] => write_buffer
    end
    context("setWriteBuffer should register a write buffer for an Edge (marginal)") do
        write_buffer = Array(Any, 0)
        setWriteBuffer(e, write_buffer)
        @fact g.write_buffers[e] => write_buffer
    end
    context("setWriteBuffer should not register a write buffer if the specified interface/edge is not in the specified graph") do
        write_buffer = Array(Any, 0)
        @fact_throws setWriteBuffer(TerminalNode().out, write_buffer)
        g2 = FactorGraph()
        @fact_throws setWriteBuffer(e, write_buffer)
    end

    # clearBuffers!
    setCurrentGraph(g)
    context("clearBuffers! should deregister all read/write buffers") do
        clearBuffers!(g)
        @fact length(g.read_buffers) => 0
        @fact length(g.write_buffers) => 0
    end
end

facts("TimeWrap integration tests") do
    # Background
    g = FactorGraph()
    node_t1 = TerminalNode()
    node_t2 = TerminalNode()
    e = Edge(node_t1, node_t2)

    # addTimeWrap
    context("addTimeWrap should not register a timewrap if the nodes do not belong to the specified graph or if the nodes are the same") do
        @fact_throws addTimeWrap(TerminalNode(), node_t1)
        @fact_throws addTimeWrap(node_t2, TerminalNode())
        @fact_throws addTimeWrap(node_t1, node_t1)
    end
    context("addTimeWrap should register a timewrap for a pair of TerminalNodes") do
        time_wraps = addTimeWrap(node_t1, node_t2)
        @fact length(time_wraps) => 1
        @fact ((node_t1, node_t2) in time_wraps) => true
    end
    context("addTimeWrap should refuse to register multiple timewraps per TerminalNode") do
        @fact_throws addTimeWrap(node_t1, node_t2)
    end

    # clearTimeWraps!
    context("clearTimeWraps! should deregister all time wraps") do
        clearTimeWraps!(g)
        @fact length(g.time_wraps) => 0
    end
end

facts("step integration tests") do
    context("step should perform a time step and handle read/write buffers") do
        # out = in + delta
        g = FactorGraph()
        node_in = TerminalNode(0.0, name="in")
        node_add = AdditionNode(name="add")
        node_delta = TerminalNode(name="delta")
        node_out = TerminalNode(name="out")
        Edge(node_in, node_add.in1, Float64)
        Edge(node_delta, node_add.in2, Float64)
        Edge(node_add.out, node_out, Float64)
        addTimeWrap(node_out, node_in)
        generateSchedule!(node_add.out)
        deltas = [1.:10.]
        setReadBuffer(node_delta, deltas)
        results = setWriteBuffer(node_add.out)
        while !isempty(deltas)
            step()
        end
        @fact results => cumsum([1.:10.])
    end
end