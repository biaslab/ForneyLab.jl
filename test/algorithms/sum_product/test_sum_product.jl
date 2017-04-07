module SumProductTest

using Base.Test
import ForneyLab: FactorGraph, Variable, GaussianMeanVariance, constant, sumProductSchedule

@testset "sumProductSchedule" begin
    FactorGraph()
    x = Variable()
    GaussianMeanVariance(x, constant(0.0), constant(1.0))

    schedule = sumProductSchedule(x)
    println(schedule)
    
    @test true == false
end

end # module