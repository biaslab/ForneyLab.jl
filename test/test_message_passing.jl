#####################
# Unit tests
#####################

facts("Message passing unit tests") do
    context("Should throw an error if one or more interfaces have no partner") do
        FixedGainNode(id=:node)
        @fact_throws execute(ForneyLab.convert(Schedule, [ForneyLab.n(:node).i[:out]]))
    end
end


#####################
# Integration tests
#####################

facts("Message passing integration tests") do
    context("execute()") do
        context("Should correctly execute a schedule and return the result of the last step") do
            initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(n(:add).i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
            setMessage!(n(:add).i[:out], Message(GaussianDistribution()))
            schedule = SumProduct.generateSchedule(n(:add).i[:in2])
            dist = ensureMVParametrization!(execute(schedule).payload)
            @fact dist --> n(:add).i[:in2].message.payload
            @fact isApproxEqual(dist.m, [2.0]) --> true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) --> true
        end

        context("Should handle post-processing of messages (sample)") do
            initializeAdditionNode(Any[Message(GaussianDistribution()), Message(GaussianDistribution()), Message(GaussianDistribution())])
            schedule = [ScheduleEntry(n(:add_node).i[:out], sumProduct!, sample)]
            @fact typeof(execute(schedule).payload) --> DeltaDistribution{Float64}
        end

        context("Should work as expeced in loopy graphs") do
            initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
            setMessage!(n(:driver).i[:out], Message(GaussianDistribution()))
            schedule = SumProduct.generateSchedule(n(:driver).i[:out])
            for count = 1:100
                execute(schedule)
            end
            @fact typeof(n(:driver).i[:out].message) --> Message{GaussianDistribution}
            @fact ensureMVParametrization!(n(:driver).i[:out].message.payload).m --> 100.0 # For stop conditions at 100 cycles deep
        end

        context("Should be called repeatedly until convergence") do
            initializeLoopyGraph(A=[1.1], B=[0.1], noise_m=0.0, noise_V=0.1)
            # Now set a breaker message and check that it works
            breaker_message = Message(GaussianDistribution(m=10.0, V=100.0))
            setMessage!(n(:driver).i[:out], breaker_message)
            prev_dist = deepcopy(breaker_message.payload)
            converged = false
            schedule = SumProduct.generateSchedule(n(:driver).i[:out])
            while !converged
                dist = ensureMVParametrization!(execute(schedule).payload)
                converged = isApproxEqual(prev_dist.m, dist.m)
                prev_dist = deepcopy(dist)
            end
            @fact isApproxEqual(n(:driver).i[:out].message.payload.m, 0.0) --> true
        end
    end

    context("clearMessage!() and clearMessages!() should clear messages") do
        initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
        setMessage!(n(:add).i[:in1], Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(n(:add).i[:out], Message(GaussianDistribution()))
        schedule = SumProduct.generateSchedule(n(:add).i[:in2])
        execute(schedule)
        clearMessage!(n(:add).i[:in2])
        @fact n(:add).i[:in2].message --> nothing
        clearMessages!(n(:add))
        @fact n(:add).i[:in1].message --> nothing
        @fact n(:add).i[:out].message --> nothing
    end
end
