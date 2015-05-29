#####################
# Integration tests
#####################

facts("Read/write buffer integration tests") do
    # setReadBuffer
    context("setReadBuffer should register a read buffer for a TerminalNode") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        read_buffer = zeros(10)
        setReadBuffer(node_t1, read_buffer)
        @fact g.read_buffers[node_t1] => read_buffer
    end
    # context("setReadBuffer should not register a read buffer if the specified node is not in the specified graph") do
    #     @fact_throws setReadBuffer(TerminalNode(), zeros(10))
    # end

    context("setReadBuffer should register a mini-batch read buffer for a TerminalNode array") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        graph = currentGraph()
        more_data = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]
        setReadBuffer(y_nodes, more_data, graph)
        @fact graph.read_buffers[y_nodes[1]] => [1.0, 2.0, 3.0]
        @fact graph.read_buffers[y_nodes[2]] => [1.0, 2.0, 3.0]
        @fact graph.read_buffers[y_nodes[3]] => [1.0, 2.0, 3.0]
    end
    
    # setWriteBuffer
    context("setWriteBuffer should register a write buffer for an Interface") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(node_t1.i[:out], write_buffer)
        @fact g.write_buffers[node_t1.i[:out]] => write_buffer
    end
    context("setWriteBuffer should register a write buffer for an Edge (marginal)") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(e, write_buffer)
        @fact g.write_buffers[e] => write_buffer
    end

    context("clearBuffers should deregister all read/write buffers") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        clearBuffers(g)
        @fact length(g.read_buffers) => 0
        @fact length(g.write_buffers) => 0
    end
end

facts("Wrap integration tests") do
    # wrap()
    context("wrap() should register a timewrap for a pair of TerminalNodes") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        wraps = wrap(node_t1, node_t2)
        @fact length(wraps) => 1
        @fact ((node_t1, node_t2) in wraps) => true
    end

    # clearWraps
    context("clearWraps should deregister all time wraps") do
        (node_t1, node_t2, e, g) = initializeBufferGraph()
        clearWraps(g)
        @fact length(g.wraps) => 0
    end
end

facts("step integration tests") do
    context("step should perform a time step and handle read/write buffers") do
        # out = in + delta
        g = FactorGraph()
        node_in = TerminalNode(DeltaDistribution(0.0), name="in")
        node_add = AdditionNode(name="add")
        node_delta = TerminalNode(name="delta")
        node_out = TerminalNode(name="out")
        Edge(node_in, node_add.i[:in1])
        Edge(node_delta, node_add.i[:in2])
        Edge(node_add.i[:out], node_out)
        wrap(node_out, node_in)
        deltas = [DeltaDistribution(n) for n in 1.:10.]
        setReadBuffer(node_delta, deltas)
        results = setWriteBuffer(node_add.i[:out])
        algo = SumProduct.Algorithm(g) # The timewraps and buffers tell the autoscheduler what should be computed
        while !isempty(deltas)
            step(algo, g)
        end
        @fact results => [DeltaDistribution(r) for r in cumsum([1.:10.])]
    end
end

facts("run() integration tests") do
    context("run() should step() until a read buffer is exhausted") do
        # out = in + delta
        g = FactorGraph()
        node_in = TerminalNode(DeltaDistribution(0.0), name="in")
        node_add = AdditionNode(name="add")
        node_delta = TerminalNode(name="delta")
        node_out = TerminalNode(name="out")
        Edge(node_in, node_add.i[:in1])
        Edge(node_delta, node_add.i[:in2])
        Edge(node_add.i[:out], node_out)
        wrap(node_out, node_in)
        deltas = [DeltaDistribution(n) for n in 1.:10.]
        setReadBuffer(node_delta, deltas)
        results = setWriteBuffer(node_add.i[:out])
        schedule = SumProduct.generateSchedule(node_add.i[:out])
        algo = Algorithm(schedule, g)
        run(algo, g)
        @fact results => [DeltaDistribution(r) for r in cumsum([1.:10.])]
    end
end