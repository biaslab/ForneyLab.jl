module WishartTest

using Test
using ForneyLab
using ForneyLab: prod!, unsafeMean, unsafeVar, unsafeDetLogMean, outboundType, isApplicable, dims, isProper, logPdf
using ForneyLab: SPWishartOutNPP, VBWishartOut
using SpecialFunctions: digamma

@testset "dims" begin
    @test dims(ProbabilityDistribution(MatrixVariate, Wishart, v=diageye(3), nu=4.0)) == (3, 3)
end

@testset "vague" begin
    @test vague(Wishart, 3) == ProbabilityDistribution(MatrixVariate, Wishart, v=huge*diageye(3), nu=3.0)
    @test vague(Wishart, (3,3)) == ProbabilityDistribution(MatrixVariate, Wishart, v=huge*diageye(3), nu=3.0)
    @test vague(Union{Gamma, Wishart}, (3,3)) == ProbabilityDistribution(MatrixVariate, Wishart, v=huge*diageye(3), nu=3.0)
    @test vague(Union{Gamma, Wishart}, ()) == ProbabilityDistribution(Univariate, Gamma, a=1.0, b=tiny)
end

@testset "isProper" begin
    @test isProper(ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=1.0)) == true
    @test isProper(ProbabilityDistribution(MatrixVariate, Wishart, v=transpose([-1.0]), nu=2.0)) == false
    @test isProper(ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=0.0)) == false
end

@testset "prod!" begin
    @test ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) * ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0),nu=2.0) == ProbabilityDistribution(MatrixVariate, Wishart, v=transpose([0.4999999999999999]),nu=2.0)
    @test ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) * ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0)) == ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0))
    @test ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0)) * ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) == ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0))
    @test_throws Exception ProbabilityDistribution(MatrixVariate, PointMass, m=transpose([-1.0])) * ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(MatrixVariate, Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == 3.0*[2.0 1.0; 1.0 2.0]
    @test unsafeVar(ProbabilityDistribution(MatrixVariate, Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == [24.0 15.0; 15.0 24.0]
    @test unsafeVar(ProbabilityDistribution(MatrixVariate, Wishart, v=diageye(3), nu=3.0)) == [6.0 3.0 3.0; 3.0 6.0 3.0; 3.0 3.0 6.0]
    @test unsafeDetLogMean(ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=1.0)) == digamma(0.5) + log(2)
    @test unsafeDetLogMean(ProbabilityDistribution(MatrixVariate, Wishart, v=eye(2), nu=2.0)) == digamma(0.5) + digamma(1) + 2*log(2)
end

@testset "log pdf" begin
    @test isapprox(logPdf(ProbabilityDistribution(MatrixVariate, Wishart, v=[3.0 1.0; 1.0 1.2], nu=6.0),[2.0 1.0; 1.0 2.0]), -8.15846321016661)
end


#-------------
# Update rules
#-------------

@testset "SPWishartOutNPP" begin
    @test SPWishartOutNPP <: SumProductRule{Wishart}
    @test outboundType(SPWishartOutNPP) == Message{Wishart}
    @test isApplicable(SPWishartOutNPP, [Nothing, Message{PointMass}, Message{PointMass}])

    @test ruleSPWishartOutNPP(nothing, Message(MatrixVariate, PointMass, m=transpose([1.0])), Message(Univariate, PointMass, m=2.0)) == Message(MatrixVariate, Wishart, v=transpose([1.0]), nu=2.0)
end

@testset "VBWishartOut" begin
    @test VBWishartOut <: NaiveVariationalRule{Wishart}
    @test outboundType(VBWishartOut) == Message{Wishart}
    @test isApplicable(VBWishartOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBWishartOut, [ProbabilityDistribution, ProbabilityDistribution, Nothing])

    @test ruleVBWishartOut(nothing, ProbabilityDistribution(MatrixVariate, PointMass, m=transpose([1.5])), ProbabilityDistribution(Univariate, PointMass, m=3.0)) == Message(MatrixVariate, Wishart, v=transpose([1.5]), nu=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)) == averageEnergy(Wishart, ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0), ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0)), ProbabilityDistribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(ProbabilityDistribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)) == differentialEntropy(ProbabilityDistribution(Univariate, Gamma, a=1.0, b=0.5))
end

end #module
