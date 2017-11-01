module WishartTest

using Base.Test
using ForneyLab
import ForneyLab: prod!, unsafeMean, unsafeVar, unsafeDetLogMean, outboundType, isApplicable, dims, isProper
import ForneyLab: VBWishartOut

@testset "dims" begin
    @test dims(MatrixVariate(Wishart, v=diageye(3), nu=4.0)) == (3, 3)
end

@testset "vague" begin
    @test vague(Wishart, 3) == MatrixVariate(Wishart, v=huge*diageye(3), nu=3.0)
end

@testset "isProper" begin
    @test isProper(MatrixVariate(Wishart, v=[1.0].', nu=1.0)) == true
    @test isProper(MatrixVariate(Wishart, v=[-1.0].', nu=2.0)) == false
    @test isProper(MatrixVariate(Wishart, v=[1.0].', nu=0.0)) == false
end

@testset "prod!" begin
    @test MatrixVariate(Wishart, v=[1.0].',nu=2.0) * MatrixVariate(Wishart, v=[1.0].',nu=2.0) == MatrixVariate(Wishart, v=[0.4999999999999999].',nu=2.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(MatrixVariate(Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == 3.0*[2.0 1.0; 1.0 2.0]
    @test unsafeVar(MatrixVariate(Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == [24.0 15.0; 15.0 24.0]
    @test unsafeVar(MatrixVariate(Wishart, v=diageye(3), nu=3.0)) == [6.0 3.0 3.0; 3.0 6.0 3.0; 3.0 3.0 6.0]
    @test unsafeDetLogMean(MatrixVariate(Wishart, v=[1.0].', nu=1.0)) == digamma(0.5) + log(2)
    @test unsafeDetLogMean(MatrixVariate(Wishart, v=eye(2), nu=2.0)) == digamma(0.5) + digamma(1) + 2*log(2)
end


#-------------
# Update rules
#-------------

@testset "VBWishartOut" begin
    @test VBWishartOut <: VariationalRule{Wishart}
    @test outboundType(VBWishartOut) == Message{Wishart}
    @test isApplicable(VBWishartOut, [Void, ProbabilityDistribution, ProbabilityDistribution]) 
    @test !isApplicable(VBWishartOut, [ProbabilityDistribution, ProbabilityDistribution, Void]) 

    @test ruleVBWishartOut(nothing, MatrixVariate(PointMass, m=[1.5].'), Univariate(PointMass, m=3.0)) == Message(MatrixVariate(Wishart, v=[1.5].', nu=3.0))
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(MatrixVariate(Wishart, v=[1.0].', nu=2.0)) == averageEnergy(Wishart, MatrixVariate(Wishart, v=[1.0].', nu=2.0), MatrixVariate(PointMass, m=[1.0].'), Univariate(PointMass, m=2.0))
    @test differentialEntropy(MatrixVariate(Wishart, v=[1.0].', nu=2.0)) == differentialEntropy(Univariate(Gamma, a=1.0, b=0.5))
end

end #module