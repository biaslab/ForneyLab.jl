#####################
# Unit tests
#####################

facts("Message passing unit tests") do
    context("calculateMessage!()") do
        context("Should throw an error if one or more interfaces have no partner") do
            node = FixedGainNode()
            scheme = InferenceScheme()
            @fact_throws calculateMessage!(node.out, scheme)
        end
    end
end


#####################
# Integration tests
#####################

facts("Message passing integration tests") do
    context("pushRequiredInbound!() should add the proper message/marginal") do
        # Composite node
        node = initializeGainEqualityCompositeNode(eye(1), true, Any[Message(DeltaDistribution(1.0)), Message(DeltaDistribution(2.0)), Message(DeltaDistribution(3.0))])
        scheme = InferenceScheme()
        @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.in1, node.out)[1], node.in1.partner.message) => true
        @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.in2, node.out)[1], node.in2.partner.message) => true

        # Not factorized
        (node, edges) = initializeGaussianNode()
        scheme = InferenceScheme()
        @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.mean, node.out)[1], node.mean.partner.message) => true
        @fact is(ForneyLab.pushRequiredInbound!(scheme, Array(Any, 0), node, node.precision, node.out)[1], node.precision.partner.message) => true
    end

    context("execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            scheme = currentScheme()
            setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
            setMessage!(add.out, Message(GaussianDistribution()))
            schedule = generateSchedule(add.in2)
            dist = ensureMVParametrization!(execute(schedule, scheme).payload)
            @fact dist => add.in2.message.payload
            @fact isApproxEqual(dist.m, [2.0]) => true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) => true
        end

        context("Should handle post-processing of messages (sample)") do
            node = initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])
            scheme = currentScheme()
            schedule = [ScheduleEntry(node.out, sumProduct!, sample)]
            @fact typeof(execute(schedule, scheme).payload) => DeltaDistribution{Array{Float64, 1}}
        end

        context("Should work as expeced in loopy graphs") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            scheme = currentScheme()
            setMessage!(driver.out, Message(GaussianDistribution()))
            schedule = generateSchedule(driver.out)
            for count = 1:100
                execute(schedule, scheme)
            end
            @fact typeof(driver.out.message) => Message{GaussianDistribution}
            @fact ensureMVParametrization!(driver.out.message.payload).m => [100.0] # For stop conditions at 100 cycles deep
        end

        context("Should be called repeatedly until convergence") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[1.1], B=[0.1], noise_m=0.0, noise_V=0.1)
            scheme = currentScheme()
            # Now set a breaker message and check that it works
            breaker_message = Message(GaussianDistribution(m=10.0, V=100.0))
            setMessage!(driver.out, breaker_message)
            prev_dist = deepcopy(breaker_message.payload)
            converged = false
            schedule = generateSchedule(driver.out)
            while !converged
                dist = ensureMVParametrization!(execute(schedule, scheme).payload)
                converged = isApproxEqual(prev_dist.m, dist.m)
                prev_dist = deepcopy(dist)
            end
            @fact isApproxEqual(driver.out.message.payload.m, [0.0]) => true
        end
    end

    context("clearMessage!() and clearMessages!() should clear messages") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
        scheme = currentScheme()
        setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(add.out, Message(GaussianDistribution()))
        schedule = generateSchedule(add.in2)
        execute(schedule, scheme)
        clearMessage!(add.in2)
        @fact add.in2.message => nothing
        clearMessages!(add)
        @fact add.in1.message => nothing
        @fact add.out.message => nothing
    end
end
