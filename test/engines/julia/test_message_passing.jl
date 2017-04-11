module JuliaMessagePassingTest

using Base.Test
import ForneyLab: FactorGraph, Variable, GaussianMeanVariance, constant, sumProductSchedule
import ForneyLab.Julia: messagePassingAlgorithm

@testset "sumProductSchedule" begin
    FactorGraph()
    x = Variable()
    GaussianMeanVariance(x, constant(0.0), constant(1.0))

    schedule = sumProductSchedule(x)
    algo = messagePassingAlgorithm(schedule)

    println(algo)

    @test true == false
end

end # module