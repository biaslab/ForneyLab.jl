#####################
# Unit tests
#####################

facts("Message passing unit tests") do
    context("Should throw an error if one or more interfaces have no partner") do
        node = FixedGainNode()
        @fact_throws execute(ForneyLab.convert(Schedule, [node.out]))
    end
end


#####################
# Integration tests
#####################

facts("Message passing integration tests") do
    context("execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
            setMessage!(add.out, Message(GaussianDistribution()))
            schedule = SumProduct.generateSchedule(add.in2)
            dist = ensureMVParametrization!(execute(schedule).payload)
            @fact dist => add.in2.message.payload
            @fact isApproxEqual(dist.m, [2.0]) => true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) => true
        end

        context("Should handle post-processing of messages (sample)") do
            node = initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])
            schedule = [ScheduleEntry(node.out, sumProduct!, sample)]
            @fact typeof(execute(schedule).payload) => DeltaDistribution{Array{Float64, 1}}
        end

        context("Should work as expeced in loopy graphs") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(driver.out, Message(GaussianDistribution()))
            schedule = SumProduct.generateSchedule(driver.out)
            for count = 1:100
                execute(schedule)
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
            schedule = SumProduct.generateSchedule(driver.out)
            while !converged
                dist = ensureMVParametrization!(execute(schedule).payload)
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
        schedule = SumProduct.generateSchedule(add.in2)
        execute(schedule)
        clearMessage!(add.in2)
        @fact add.in2.message => nothing
        clearMessages!(add)
        @fact add.in1.message => nothing
        @fact add.out.message => nothing
    end
end
