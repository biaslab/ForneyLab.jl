#####################
# Unit tests
#####################

facts("Message passing unit tests") do
    context("calculateMessage!()") do
        context("Should throw an error if the specified interface does not belong to the specified node") do
            (node1, node2) = initializePairOfNodes()
            @fact_throws calculateMessage!(node1.out, node2)
        end

        context("Should throw an error if one or more interfaces have no partner") do
            node = FixedGainNode()
            @fact_throws calculateMessage!(node.out)
        end
    end
end


#####################
# Integration tests
#####################

facts("Message passing integration tests") do
    context("pushRequiredInbound!() should add the proper message/marginal") do
        # Composite node
        node = initializeGainEqualityCompositeNode(eye(1), true, Any[Message(1.0), Message(2.0), Message(3.0)])
        graph = getCurrentGraph()
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.in1, node.out)[1], node.in1.partner.message) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.in2, node.out)[1], node.in2.partner.message) => true

        # Not factorized
        (node, edges) = initializeGaussianNode()
        graph = getCurrentGraph()
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.mean, node.out)[1], node.mean.partner.message) => true 
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.precision, node.out)[1], node.precision.partner.message) => true

        # Mean field factorized Gaussian node
        (node, edges) = initializeGaussianNode()
        graph = getCurrentGraph()
        factorizeMeanField!(graph)
        setUninformativeMarginals!(graph)
        sg_mean = getSubgraph(node.mean.edge)
        sg_prec = getSubgraph(node.precision.edge)
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.mean, node.out)[1], graph.approximate_marginals[(node, sg_mean)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.precision, node.out)[1], graph.approximate_marginals[(node, sg_prec)]) => true

        # Mean field factorized linear node
        node = initializeLinearCompositeNode()
        graph = getCurrentGraph()
        factorizeMeanField!(graph)
        setUninformativeMarginals!(graph)
        sg_a = getSubgraph(node.slope.edge)
        sg_b = getSubgraph(node.offset.edge)
        sg_gam = getSubgraph(node.noise.edge)
        sg_x = getSubgraph(node.in1.edge)
        sg_y = getSubgraph(node.out.edge)
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.slope, node.out)[1], graph.approximate_marginals[(node, sg_a)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.offset, node.out)[1], graph.approximate_marginals[(node, sg_b)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.noise, node.out)[1], graph.approximate_marginals[(node, sg_gam)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.in1, node.out)[1], graph.approximate_marginals[(node, sg_x)]) => true

        # Structurally factorized
        (node, edges) = initializeGaussianNode()
        graph = getCurrentGraph()
        factorize!(node.out.edge)
        setUninformativeMarginals!(graph)
        sg_mean_prec = getSubgraph(node.mean.edge)
        sg_out = getSubgraph(node.out.edge)
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.mean, node.out)[1], graph.approximate_marginals[(node, sg_mean_prec)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.precision, node.out)[1], graph.approximate_marginals[(node, sg_mean_prec)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.precision, node.mean)[1], node.precision.partner.message) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.mean, node.precision)[1], node.mean.partner.message) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.out, node.mean)[1], graph.approximate_marginals[(node, sg_out)]) => true
        @fact is(ForneyLab.pushRequiredInbound!(graph, Array(Any, 0), node, node.out, node.precision)[1], graph.approximate_marginals[(node, sg_out)]) => true
    end

    context("calculateMessage!()") do
        context("Should return and write back an output message") do
            (gain, terminal) = initializePairOfNodes(A=[2.0], msg_gain_1=nothing, msg_gain_2=nothing, msg_terminal=Message(3.0))
            Edge(terminal.out, gain.in1, Float64, Array{Float64, 2})
            Edge(gain.out, MockNode().out, Array{Float64, 2})
            gain.out.message_payload_type = Array{Float64, 2} # Expect a matrix
            @fact gain.out.message => nothing
            # Request message on node for which the input is unknown
            msg = calculateMessage!(gain.out)
            @fact msg => gain.out.message # Returned message should be identical to message stored on interface.
            @fact typeof(gain.out.message.payload) => Array{Float64, 2}
            @fact gain.out.message.payload => reshape([6.0], 1, 1)
        end

        context("Should recursively calculate required inbound message") do
            # Define three nodes in series
            (node1, node2, node3) = initializeChainOfNodes()
            @fact node3.out.message => nothing
            # Request message on node for which the input is unknown
            node3.out.message_payload_type = Array{Float64, 2} # Expect a matrix
            calculateMessage!(node3.out)
            @fact typeof(node3.out.message.payload) => Array{Float64, 2}
            @fact node3.out.message.payload => reshape([12.0], 1, 1)
        end

        context("Should throw an error when there is an unbroken loop") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            @fact_throws calculateMessage!(driver.out)
        end
    end

    context("executeSchedule()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
            setMessage!(add.out, Message(GaussianDistribution()))
            schedule = generateSchedule(add.in2)
            dist = ensureMVParametrization!(executeSchedule(schedule).payload)
            @fact dist => add.in2.message.payload
            @fact isApproxEqual(dist.m, [2.0]) => true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) => true
        end

        context("Should work as expeced in loopy graphs") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(driver.out, Message(GaussianDistribution()))
            schedule = generateSchedule(driver.out)
            for count = 1:100
                executeSchedule(schedule)
            end
            @fact typeof(driver.out.message) => Message{GaussianDistribution}
            @fact ensureMVParametrization!(driver.out.message.payload).m => [100.0] # For stop conditions at 100 cycles deep
        end

        context("Should be called repeatedly until convergence") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[1.1], B=[0.1], noise_m=0.0, noise_V=0.1)
            # Now set a breaker message and check that it works
            breaker_message = Message(GaussianDistribution(m=10.0, V=100.0))
            setMessage!(driver.out, breaker_message)
            prev_dist = deepcopy(breaker_message.payload)
            converged = false
            schedule = generateSchedule(driver.out)
            while !converged
                dist = ensureMVParametrization!(executeSchedule(schedule).payload)
                converged = isApproxEqual(prev_dist.m, dist.m)
                prev_dist = deepcopy(dist)
            end
            @fact isApproxEqual(driver.out.message.payload.m, [0.0]) => true
        end
    end

    context("clearMessage!() and clearMessages!() should clear messages") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
        setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(add.out, Message(GaussianDistribution()))
        schedule = generateSchedule(add.in2)
        executeSchedule(schedule)
        clearMessage!(add.in2)
        @fact add.in2.message => nothing
        clearMessages!(add)
        @fact add.in1.message => nothing
        @fact add.out.message => nothing
    end
end