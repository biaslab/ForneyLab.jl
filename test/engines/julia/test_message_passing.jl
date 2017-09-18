module JuliaMessagePassingTest

using Base.Test
using ForneyLab

@testset "Julia.messagePassingAlgorithm" begin
    @testset "Read from data" begin
        g = FactorGraph()
        x = Variable()
        m = Variable()
        GaussianMeanVariance(m, constant(0.0), constant(1.0))
        GaussianMeanVariance(x, m, constant(1.0))
        placeholder(x, :x)

        schedule = sumProductSchedule(m)
        algo = ForneyLab.Julia.messagePassingAlgorithm(schedule, m)

        @test contains(algo, "messages = Array{Message}(2)")
        @test contains(algo, "messages[1] = rule(ForneyLab.SPGaussianMeanVariancePPV, [Message{PointMass}(0.0), Message{PointMass}(1.0), nothing])")
        @test contains(algo, "messages[2] = rule(ForneyLab.SPGaussianMeanVarianceVPP, [nothing, Message{PointMass}(1.0), Message{PointMass}(data[:x])])")
        @test contains(algo, "marginals[:variable_2] = messages[1].dist * messages[2].dist")
    end

    @testset "Model with dangling edge" begin
        FactorGraph()
        x = Variable()
        GaussianMeanVariance(x, constant(0.0), constant(1.0))

        schedule = sumProductSchedule(x)
        algo = ForneyLab.Julia.messagePassingAlgorithm(schedule, x)

        @test contains(algo, "messages = Array{Message}(1)")
        @test contains(algo, "messages[1] = rule(ForneyLab.SPGaussianMeanVariancePPV, [Message{PointMass}(0.0), Message{PointMass}(1.0), nothing])")
        @test contains(algo, "marginals[:variable_1] = messages[1].dist")
    end
end

end # module