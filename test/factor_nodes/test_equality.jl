module EqualityTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPEqualityGaussian


#-------------
# Update rules
#-------------

@testset "SPEqualityGaussian" begin
    @test SPEqualityGaussian <: SumProductRule{Equality}
    @test outboundType(SPEqualityGaussian) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPEqualityGaussian, [Message{GaussianWeightedMeanPrecision}, Message{GaussianWeightedMeanPrecision}, Void]) 
    @test isApplicable(SPEqualityGaussian, [Void, Message{GaussianWeightedMeanPrecision}, Message{GaussianMeanVariance}]) 
end

# TODO: Add more tests
@testset "ruleSPEqualityGaussian" begin
    @test ruleSPEqualityGaussian(Message(GaussianWeightedMeanPrecision, xi=1.0, w=2.0), Message(GaussianWeightedMeanPrecision, xi=3.0, w=4.0), nothing) == Message(GaussianWeightedMeanPrecision, xi=4.0, w=6.0)
end

end #module