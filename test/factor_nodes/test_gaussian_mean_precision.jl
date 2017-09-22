module GaussianMeanPrecisionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, VBGaussianMeanPrecision1, VBGaussianMeanPrecision2, VBGaussianMeanPrecision3


#-------------
# Update rules
#-------------

@testset "VBGaussianMeanPrecision1" begin
    @test VBGaussianMeanPrecision1 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision1) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecision1, GaussianMeanPrecision, 1) 
    @test !isApplicable(VBGaussianMeanPrecision1, GaussianMeanPrecision, 2) 

    @test ruleVBGaussianMeanPrecision1(nothing, ProbabilityDistribution(Gamma, a=1.0, b=2.0), ProbabilityDistribution(Gaussian, m=3.0, v=4.0)) == Message(Gaussian, m=3.0, w=0.5)
end

@testset "VBGaussianMeanPrecision2" begin
    @test VBGaussianMeanPrecision2 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision2) == Message{Gamma}
    @test isApplicable(VBGaussianMeanPrecision2, GaussianMeanPrecision, 2) 

    @test ruleVBGaussianMeanPrecision2(ProbabilityDistribution(Gaussian, m=1.0, v=2.0), Void, ProbabilityDistribution(Gaussian, m=3.0, v=4.0)) == Message(Gamma, a=1.5, b=0.5*(2.0 + 4.0 + (3.0 - 1.0)^2))
end

@testset "VBGaussianMeanPrecision3" begin
    @test VBGaussianMeanPrecision3 <: VariationalRule{GaussianMeanPrecision}
    @test outboundType(VBGaussianMeanPrecision3) == Message{Gaussian}
    @test isApplicable(VBGaussianMeanPrecision3, GaussianMeanPrecision, 3) 

    @test ruleVBGaussianMeanPrecision3(ProbabilityDistribution(Gaussian, m=3.0, v=4.0), ProbabilityDistribution(Gamma, a=1.0, b=2.0), nothing) == Message(Gaussian, m=3.0, w=0.5)
end

end #module