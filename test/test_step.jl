#####################
# Integration tests
#####################

facts("Read/write buffer integration tests") do
    # setReadBuffer
    context("setReadBuffer should register a read buffer for a TerminalNode") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        read_buffer = zeros(10)
        setReadBuffer(node_t1, read_buffer)
        @fact s.read_buffers[node_t1] => read_buffer
    end
    # context("setReadBuffer should not register a read buffer if the specified node is not in the specified graph") do
    #     @fact_throws setReadBuffer(TerminalNode(), zeros(10))
    # end

    context("setReadBuffer should register a mini-batch read buffer for a TerminalNode array") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        scheme = currentScheme()
        more_data = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]
        setReadBuffer(y_nodes, more_data, scheme)
        @fact scheme.read_buffers[y_nodes[1]] => [1.0, 2.0, 3.0]
        @fact scheme.read_buffers[y_nodes[2]] => [1.0, 2.0, 3.0]
        @fact scheme.read_buffers[y_nodes[3]] => [1.0, 2.0, 3.0]
    end
    
    # setWriteBuffer
    context("setWriteBuffer should register a write buffer for an Interface") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(node_t1.out, write_buffer)
        @fact s.write_buffers[node_t1.out] => write_buffer
    end
    context("setWriteBuffer should register a write buffer for an Edge (marginal)") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(e, write_buffer)
        @fact s.write_buffers[e] => write_buffer
    end

    context("clearBuffers should deregister all read/write buffers") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        clearBuffers(s)
        @fact length(s.read_buffers) => 0
        @fact length(s.write_buffers) => 0
    end
end

facts("TimeWrap integration tests") do
    # setTimeWrap
    context("setTimeWrap should register a timewrap for a pair of TerminalNodes") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        time_wraps = setTimeWrap(node_t1, node_t2)
        @fact length(time_wraps) => 1
        @fact ((node_t1, node_t2) in time_wraps) => true
    end

    # clearTimeWraps
    context("clearTimeWraps should deregister all time wraps") do
        (node_t1, node_t2, e, s) = initializeBufferGraph()
        clearTimeWraps(s)
        @fact length(s.time_wraps) => 0
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
        Edge(node_in, node_add.in1)
        Edge(node_delta, node_add.in2)
        Edge(node_add.out, node_out)
        s = InferenceScheme()
        setTimeWrap(node_out, node_in)
        generateSchedule!(node_add.out)
        deltas = [DeltaDistribution(n) for n in 1.:10.]
        setReadBuffer(node_delta, deltas)
        results = setWriteBuffer(node_add.out)
        while !isempty(deltas)
            step(s)
        end
        @fact results => [DeltaDistribution(r) for r in cumsum([1.:10.])]
    end

    context("step should accept and execute a number of iterations for VMP") do
        data = [GaussianDistribution(m=2.0, V=tiny())]
        g = FactorGraph()
        g_node = GaussianNode(form="precision")
        t_out = TerminalNode(name="t_out")
        t_mean = TerminalNode(GaussianDistribution(), name="t_mean")
        t_var = TerminalNode(GammaDistribution(), name="t_var")
        Edge(g_node.out, t_out, GaussianDistribution)
        Edge(t_mean, g_node.mean, GaussianDistribution)
        Edge(t_var, g_node.precision, GammaDistribution)
        s = InferenceScheme()

        setReadBuffer(t_out, data)
        mean_out = setWriteBuffer(g_node.mean.edge)
        prec_out = setWriteBuffer(g_node.precision.edge)
        factorize!(s)
        for subgraph in s.factorization
            generateSchedule!(subgraph)
        end
        setVagueQDistributions!()
        step(s, n_iterations=10)
        @fact round(mean_out[end].W[1,1], 2) => 1.79
        @fact round(mean_out[end].xi[1], 2) => 1.57
        @fact prec_out[end].a => 1.5
        @fact round(prec_out[end].b, 2) => 1.91
    end
end