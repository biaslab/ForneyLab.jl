#####################
# Integration tests
#####################

facts("Read/write buffer integration tests") do
    # setReadBuffer
    context("setReadBuffer should register a read buffer for a TerminalNode") do
        g = initializeBufferGraph()
        read_buffer = zeros(10)
        setReadBuffer(n(:node_t1), read_buffer)
        @fact g.read_buffers[n(:node_t1)] => read_buffer
    end
    # context("setReadBuffer should not register a read buffer if the specified node is not in the specified graph") do
    #     @fact_throws setReadBuffer(TerminalNode(), zeros(10))
    # end

    context("setReadBuffer should register a mini-batch read buffer for a TerminalNode array") do
        data = [1.0, 1.0, 1.0]
        initializeGaussianNodeChain(data)
        graph = currentGraph()
        more_data = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]
        setReadBuffer([n(:y1), n(:y2), n(:y3)], more_data, graph)
        @fact graph.read_buffers[n(:y1)] => [1.0, 2.0, 3.0]
        @fact graph.read_buffers[n(:y2)] => [1.0, 2.0, 3.0]
        @fact graph.read_buffers[n(:y3)] => [1.0, 2.0, 3.0]
    end
    
    # setWriteBuffer
    context("setWriteBuffer should register a write buffer for an Interface") do
        g = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(n(:node_t1).i[:out], write_buffer)
        @fact g.write_buffers[n(:node_t1).i[:out]] => write_buffer
    end
    context("setWriteBuffer should register a write buffer for an Edge (marginal)") do
        g = initializeBufferGraph()
        write_buffer = Array(Any, 0)
        setWriteBuffer(e(:e), write_buffer)
        @fact g.write_buffers[e(:e)] => write_buffer
    end

    context("clearBuffers should deregister all read/write buffers") do
        g = initializeBufferGraph()
        clearBuffers(g)
        @fact length(g.read_buffers) => 0
        @fact length(g.write_buffers) => 0
    end
end

facts("Wrap integration tests") do
    # Wrap()
    context("Wrap() should register a timewrap for a pair of TerminalNodes") do
        g = initializeBufferGraph()
        wraps = Wrap(n(:node_t1), n(:node_t2))
        @fact length(wraps) => 1
        @fact ((n(:node_t1), n(:node_t2)) in wraps) => true
    end

    # clearWraps
    context("clearWraps should deregister all time wraps") do
        g = initializeBufferGraph()
        clearWraps(g)
        @fact length(g.wraps) => 0
    end
end

facts("step integration tests") do
    context("step should perform a time step and handle read/write buffers") do
        # out = in + delta
        g = FactorGraph()
        TerminalNode(DeltaDistribution(0.0), id=:in)
        AdditionNode(id=:add)
        TerminalNode(id=:delta)
        TerminalNode(id=:out)
        Edge(n(:in), n(:add).i[:in1])
        Edge(n(:delta), n(:add).i[:in2])
        Edge(n(:add).i[:out], n(:out))
        Wrap(n(:out), n(:in))
        deltas = [DeltaDistribution(n) for n in 1.:10.]
        setReadBuffer(n(:delta), deltas)
        results = setWriteBuffer(n(:add).i[:out])
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
        TerminalNode(DeltaDistribution(0.0), id=:in)
        AdditionNode(id=:add)
        TerminalNode(id=:delta)
        TerminalNode(id=:out)
        Edge(n(:in), n(:add).i[:in1])
        Edge(n(:delta), n(:add).i[:in2])
        Edge(n(:add).i[:out], n(:out))
        Wrap(n(:out), n(:in))
        deltas = [DeltaDistribution(n) for n in 1.:10.]
        setReadBuffer(n(:delta), deltas)
        results = setWriteBuffer(n(:add).i[:out])
        schedule = SumProduct.generateSchedule(n(:add).i[:out])
        algo = Algorithm(schedule, g)
        run(algo, g)
        @fact results => [DeltaDistribution(r) for r in cumsum([1.:10.])]
    end
end