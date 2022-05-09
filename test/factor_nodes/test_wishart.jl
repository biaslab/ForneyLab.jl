module WishartTest

using Test
using ForneyLab
using ForneyLab: prod!, unsafeMean, unsafeVar, unsafeDetLogMean, outboundType, isApplicable, dims, isProper, logPdf, naturalParams, standardDistribution
using ForneyLab: SPWishartOutNPP, VBWishartOut
using SpecialFunctions: digamma

@testset "dims" begin
    @test dims(Distribution(MatrixVariate, Wishart, v=diageye(3), nu=4.0)) == (3,3)
end

@testset "vague" begin
    @test vague(Wishart, (3,3)) == Distribution(MatrixVariate, Wishart, v=huge*diageye(3), nu=3.0)
    @test vague(Union{Gamma, Wishart}, (3,3)) == Distribution(MatrixVariate, Wishart, v=huge*diageye(3), nu=3.0)
    @test vague(Union{Gamma, Wishart}, ()) == Distribution(Univariate, Gamma, a=1.0, b=tiny)
end

@testset "isProper" begin
    @test isProper(Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=1.0)) == true
    @test isProper(Distribution(MatrixVariate, Wishart, v=transpose([-1.0]), nu=2.0)) == false
    @test isProper(Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=0.0)) == false
end

@testset "prod!" begin
    @test Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) * Distribution(MatrixVariate, Wishart, v=mat(1.0),nu=2.0) == Distribution(MatrixVariate, Wishart, v=transpose([0.4999999999999999]),nu=2.0)
    @test Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) * Distribution(MatrixVariate, PointMass, m=mat(1.0)) == Distribution(MatrixVariate, PointMass, m=mat(1.0))
    @test Distribution(MatrixVariate, PointMass, m=mat(1.0)) * Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0) == Distribution(MatrixVariate, PointMass, m=mat(1.0))
    @test_throws Exception Distribution(MatrixVariate, PointMass, m=transpose([-1.0])) * Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(MatrixVariate, Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == 3.0*[2.0 1.0; 1.0 2.0]
    @test unsafeVar(Distribution(MatrixVariate, Wishart, v=[2.0 1.0; 1.0 2.0], nu=3.0)) == [24.0 15.0; 15.0 24.0]
    @test unsafeVar(Distribution(MatrixVariate, Wishart, v=diageye(3), nu=3.0)) == [6.0 3.0 3.0; 3.0 6.0 3.0; 3.0 3.0 6.0]
    @test unsafeDetLogMean(Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=1.0)) == digamma(0.5) + log(2)
    @test unsafeDetLogMean(Distribution(MatrixVariate, Wishart, v=eye(2), nu=2.0)) == digamma(0.5) + digamma(1) + 2*log(2)
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(MatrixVariate, Wishart, v=[3.0 1.0; 1.0 1.2], nu=6.0),[2.0 1.0; 1.0 2.0]), -8.15846321016661)
end

@testset "natural parameters" begin
    d = Distribution(MatrixVariate, Wishart, v=[1.0 0.5; 0.5 2.0], nu=4.0)
    η = naturalParams(d)
    s = standardDistribution(MatrixVariate, Wishart, η=η)
    @test isapprox(d.params[:v], s.params[:v]) # Test conversion consistency
    @test d.params[:nu] == s.params[:nu]

    x = [[3.0 0.1; 0.1 2.0], [1.0 0.3; 0.3 4.0]]
    d_x = logPdf.([d], x)
    η_x = logPdf.(MatrixVariate, Wishart, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
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
    @test isApplicable(VBWishartOut, [Nothing, Distribution, Distribution])
    @test !isApplicable(VBWishartOut, [Distribution, Distribution, Nothing])

    @test ruleVBWishartOut(nothing, Distribution(MatrixVariate, PointMass, m=transpose([1.5])), Distribution(Univariate, PointMass, m=3.0)) == Message(MatrixVariate, Wishart, v=transpose([1.5]), nu=3.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test differentialEntropy(Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)) == averageEnergy(Wishart, Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0), Distribution(MatrixVariate, PointMass, m=mat(1.0)), Distribution(Univariate, PointMass, m=2.0))
    @test differentialEntropy(Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=2.0)) == differentialEntropy(Distribution(Univariate, Gamma, a=1.0, b=0.5))
end

end #module
