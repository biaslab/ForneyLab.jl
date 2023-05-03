module PoissonTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, unsafeMean, unsafeVar, slug, isProper, naturalParams, standardDistribution
using ForneyLab: VBPoissonOut, VBPoissonL, SPPoissonOutNP, SPPoissonLPN

@testset "Poisson Distribution construction" begin
    @test Distribution(Univariate, Poisson, l=2.0) == Distribution{Univariate, Poisson}(Dict(:l=>2.0))
    @test Distribution(Poisson, l=2.0) == Distribution{Univariate, Poisson}(Dict(:l=>2.0))
    @test Distribution(Poisson) == Distribution{Univariate, Poisson}(Dict(:l=>1.0))
    @test_throws MethodError Distribution(Multivariate, Poisson, l=2.0)
end

@testset "Poisson Message construction" begin
    @test Message(Univariate, Poisson, l=2.0) == Message{Poisson, Univariate}(Distribution{Univariate, Poisson}(Dict(:l=>2.0))    )
    @test Message(Poisson, l=2.0) == Message{Poisson, Univariate}(Distribution{Univariate, Poisson}(Dict(:l=>2.0))    )
    @test Message(Poisson) == Message{Poisson, Univariate}(Distribution{Univariate, Poisson}(Dict(:l=>1.0)))
    @test_throws MethodError Message(Multivariate, Poisson, l=2.0)
end

@testset "dims" begin
    @test dims(Distribution(Poisson)) == ()
end

@testset "slug" begin
    @test slug(Poisson) == "Poisson"
end

@testset "vague" begin
    @test vague(Poisson) == Distribution(Poisson, l=huge)
end

@testset "isProper" begin
    @test isProper(Distribution(Poisson, l=1.0))
    @test !isProper(Distribution(Poisson, l=0.0))
    @test !isProper(Distribution(Poisson, l=-1.0))
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(Poisson, l=2.0)) == 2.0
    @test unsafeVar(Distribution(Poisson, l=2.0)) == 2.0
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Poisson, l=2.5),1), -1.583709268125845)
end

@testset "natural parameters" begin
    d = Distribution(Univariate, Poisson, l=2.0)
    η = naturalParams(d)
    s = standardDistribution(Univariate, Poisson, η=η)
    @test d.params[:l] == s.params[:l] # Test conversion consistency

    x = [1, 4, 8]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Univariate, Poisson, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end


#-------------
# Update rules
#-------------

@testset "SPPoissonOutNP" begin
    @test SPPoissonOutNP <: SumProductRule{Poisson}
    @test outboundType(SPPoissonOutNP) == Message{Poisson}
    @test isApplicable(SPPoissonOutNP, [Nothing, Message{PointMass}])

    @test ruleSPPoissonOutNP(nothing, Message(Univariate, PointMass, m=2.0)) == Message(Poisson, l=2.0)
end

@testset "SPPoissonLPN" begin
    @test SPPoissonLPN <: SumProductRule{Poisson}
    @test outboundType(SPPoissonLPN) == Message{Gamma}
    @test isApplicable(SPPoissonLPN, [Message{PointMass}, Nothing])

    @test ruleSPPoissonLPN(Message(Univariate, PointMass, m=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
end

@testset "VBPoissonOut" begin
    @test VBPoissonOut <: NaiveVariationalRule{Poisson}
    @test outboundType(VBPoissonOut) == Message{Poisson}
    @test isApplicable(VBPoissonOut, [Nothing, Distribution])

    @test ruleVBPoissonOut(nothing, Distribution(Gamma, a=2.0, b=3.0)) == Message(Poisson, l=0.5087350371986222)
end

@testset "VBPoissonL" begin
    @test VBPoissonL <: NaiveVariationalRule{Poisson}
    @test outboundType(VBPoissonL) == Message{Gamma}
    @test isApplicable(VBPoissonL, [Distribution, Nothing])

    @test ruleVBPoissonL(Distribution(Poisson, l=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
    @test ruleVBPoissonL(Distribution(Univariate, PointMass, m=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
end

@testset "averageEnergy and differentialEntropy" begin
    @test isapprox(differentialEntropy(Distribution(Poisson, l=1.0)), averageEnergy(Poisson, Distribution(Poisson, l=1.0), Distribution(Univariate, PointMass, m=1.0)))
    @test isapprox(differentialEntropy(Distribution(Poisson, l=10.0)), averageEnergy(Poisson, Distribution(Poisson, l=10.0), Distribution(Univariate, PointMass, m=10.0)))
    @test isapprox(differentialEntropy(Distribution(Poisson, l=100.0)), averageEnergy(Poisson, Distribution(Poisson, l=100.0), Distribution(Univariate, PointMass, m=100.0)), atol=0.1)

    @test averageEnergy(Poisson, Distribution(PointMass, m=1.0), Distribution(Univariate, PointMass, m=1.0)) == 1.0
    @test averageEnergy(Poisson, Distribution(PointMass, m=2.0), Distribution(Univariate, PointMass, m=1.0)) == 1.0 + log(2)
end

end # module
