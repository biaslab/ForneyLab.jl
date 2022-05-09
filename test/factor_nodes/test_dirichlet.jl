module DirichletTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable, prod!, unsafeMean, unsafeLogMean, unsafeVar, vague, dims, naturalParams, standardDistribution
using ForneyLab: SPDirichletOutNP, VBDirichletOut, VBDirichletIn1
using SpecialFunctions: digamma

@testset "Dirichlet Distribution and Message construction" begin
    @test Distribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == Distribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test Distribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]) == Distribution{MatrixVariate, Dirichlet}(Dict(:a=>[2.0 3.0; 4.0 5.0]))
    @test_throws Exception Distribution(Univariate, Dirichlet)
    @test Distribution(Dirichlet, a=[2.0, 3.0]) == Distribution{Multivariate, Dirichlet}(Dict(:a=>[2.0, 3.0]))
    @test Distribution(Dirichlet) == Distribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0]))
    @test Message(Dirichlet) == Message{Dirichlet, Multivariate}(Distribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test Message(Multivariate, Dirichlet) == Message{Dirichlet, Multivariate}(Distribution{Multivariate, Dirichlet}(Dict(:a=>[1.0, 1.0, 1.0])))
    @test Message(MatrixVariate, Dirichlet) == Message{Dirichlet, MatrixVariate}(Distribution{MatrixVariate, Dirichlet}(Dict(:a=>ones(3,3))))
    @test_throws Exception Message(Univariate, Dirichlet)
end

@testset "dims" begin
    @test dims(Distribution(Multivariate, Dirichlet, a=[2.0, 2.0, 2.0])) == (3,)
    @test dims(Distribution(MatrixVariate, Dirichlet, a=[2.0 2.0 2.0; 2.0 2.0 2.0])) == (2,3)
end

@testset "vague" begin
    @test vague(Dirichlet, (3,)) == Distribution(Dirichlet, a=ones(3))
    @test vague(Dirichlet, (2,3)) == Distribution(MatrixVariate, Dirichlet, a=ones(2,3))
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(Distribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.5, 0.5]
    @test unsafeMean(Distribution(MatrixVariate, Dirichlet, a=[2.0 1.0; 2.0 3.0])) == [0.5 0.25; 0.5 0.75]
    @test unsafeLogMean(Distribution(Multivariate, Dirichlet, a=[2.0, 3.0])) == [digamma(2.0), digamma(3.0)] .- digamma(5.0)
    @test unsafeLogMean(Distribution(MatrixVariate, Dirichlet, a=[2.0 4.0; 3.0 5.0])) == [digamma(2.0) digamma(4.0); digamma(3.0) digamma(5.0)] - [digamma(5.0) digamma(9.0); digamma(5.0) digamma(9.0)]
    @test unsafeVar(Distribution(Multivariate, Dirichlet, a=[2.0, 2.0])) == [0.05, 0.05]
end

@testset "prod!" begin
    # Multivariate
    @test Distribution(Multivariate, Dirichlet, a=[2.0, 2.0]) * Distribution(Multivariate, Dirichlet, a=[2.0, 3.0]) == Distribution(Multivariate, Dirichlet, a=[3.0, 4.0])
    @test Distribution(Multivariate, Dirichlet, a=[2.0, 2.0]) * Distribution(Univariate, Beta, a=2.0, b=3.0) == Distribution(Multivariate, Dirichlet, a=[3.0, 4.0])
    @test Distribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) * Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) == Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])
    @test Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1]) * Distribution(Multivariate, Dirichlet, a=[1.0, 2.0, 3.0]) == Distribution(Multivariate, PointMass, m=[0.1, 0.8, 0.1])

    # MatrixVariate
    @test Distribution(MatrixVariate, Dirichlet, a=[2.0 2.0; 3.0 3.0]) * Distribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]) == Distribution(MatrixVariate, Dirichlet, a=[3.0 4.0; 6.0 7.0])
    @test Distribution(MatrixVariate, Dirichlet, a=[1.0 3.0; 2.0 4.0]) * Distribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7]) == Distribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7])
    @test Distribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7]) * Distribution(MatrixVariate, Dirichlet, a=[1.0 3.0; 2.0 4.0]) == Distribution(MatrixVariate, PointMass, m=[0.1 0.3; 0.9 0.7])
end

@testset "log pdf" begin
    @test isapprox(logPdf(Distribution(Multivariate, Dirichlet, a=[0.2,3.0,1.5]),[2,3,7]), 3.2556382883760024)
    @test isapprox(logPdf(Distribution(MatrixVariate, Dirichlet, a=[0.2 1.4; 3.0 1.8]),[2 7; 3 3]), 3.0442561618507087)
end

@testset "natural parameters" begin
    # Multivariate
    d = Distribution(Multivariate, Dirichlet, a=[1.5, 3.0, 5.0])
    η = naturalParams(d)
    s = standardDistribution(Multivariate, Dirichlet, η=η)
    @test d.params[:a] == s.params[:a] # Test conversion consistency

    x = [[0.2, 0.6, 0.2], [0.1, 0.4, 0.5]]
    d_x = logPdf.([d], x)
    η_x = logPdf.(Multivariate, Dirichlet, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency

    # MatrixVariate
    d = Distribution(MatrixVariate, Dirichlet, a=[1.5 3.0; 5.0 4.0; 4.0 6.0])
    η = naturalParams(d)
    s = standardDistribution(MatrixVariate, Dirichlet, η=η, dims=(3,2))
    @test d.params[:a] == s.params[:a] # Test conversion consistency

    x = [[0.2 0.1; 0.6 0.4; 0.2 0.5], [0.1 0.4; 0.5 0.1; 0.4 0.5]]
    d_x = logPdf.([d], x)
    η_x = logPdf.(MatrixVariate, Dirichlet, x; η=η)
    @test isapprox(d_x, η_x) # Test pdf consistency
end


#-------------
# Update rules
#-------------

@testset "SPDirichletOutNP" begin
    @test SPDirichletOutNP <: SumProductRule{Dirichlet}
    @test outboundType(SPDirichletOutNP) == Message{Dirichlet}
    @test isApplicable(SPDirichletOutNP, [Nothing, Message{PointMass}])

    @test ruleSPDirichletOutNP(nothing, Message(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
    @test ruleSPDirichletOutNP(nothing, Message(MatrixVariate, PointMass, m=[2.0 3.0; 4.0 5.0])) == Message(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])
end

@testset "VBDirichletOut" begin
    @test VBDirichletOut <: NaiveVariationalRule{Dirichlet}
    @test outboundType(VBDirichletOut) == Message{Dirichlet}
    @test isApplicable(VBDirichletOut, [Nothing, Distribution])

    @test ruleVBDirichletOut(nothing, Distribution(Multivariate, PointMass, m=[2.0, 3.0])) == Message(Multivariate, Dirichlet, a=[2.0, 3.0])
    @test ruleVBDirichletOut(nothing, Distribution(MatrixVariate, PointMass, m=[2.0 3.0; 4.0 5.0])) == Message(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])
end

@testset "VBDirichletIn1" begin
    @test VBDirichletIn1 <: NaiveVariationalRule{Dirichlet}
    @test outboundType(VBDirichletIn1) == Message{Function}
    @test isApplicable(VBDirichletIn1, [Distribution, Nothing])
end

@testset "averageEnergy and differentialEntropy" begin
    # Multivariate
    @test isapprox(differentialEntropy(Distribution(Multivariate, Dirichlet, a=[2.0, 3.0])), averageEnergy(Dirichlet, Distribution(Multivariate, Dirichlet, a=[2.0, 3.0]), Distribution(Multivariate, PointMass, m=[2.0, 3.0])))
    @test isapprox(averageEnergy(Dirichlet, Distribution(Multivariate, Dirichlet, a=[4.0, 5.0]), Distribution(Multivariate, PointMass, m=[2.0, 3.0])), averageEnergy(Beta, Distribution(Univariate, Beta, a=4.0, b=5.0), Distribution(Univariate, PointMass, m=2.0), Distribution(Univariate, PointMass, m=3.0)))
    @test isapprox(averageEnergy(Dirichlet,Distribution(Multivariate,Dirichlet,a=[2.0,3.0,4.0]),Distribution(Multivariate,SampleList,s=[[1.0,2.0,1.0],[3.,3.,1.],[2.,2.,2]],w=[0.1,0.7,0.3])),0.5093876870003795)

    # MatrixVariate
    @test differentialEntropy(Distribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0])) == differentialEntropy(Distribution(Multivariate, Dirichlet, a=[2.0, 4.0])) + differentialEntropy(Distribution(Multivariate, Dirichlet, a=[3.0, 5.0]))
    @test averageEnergy(Dirichlet, Distribution(MatrixVariate, Dirichlet, a=[2.0 3.0; 4.0 5.0]), Distribution(MatrixVariate, PointMass, m=[6.0 7.0; 8.0 9.0])) == averageEnergy(Dirichlet, Distribution(Multivariate, Dirichlet, a=[2.0, 4.0]), Distribution(Multivariate, PointMass, m=[6.0, 8.0])) + averageEnergy(Dirichlet, Distribution(Multivariate, Dirichlet, a=[3.0, 5.0]), Distribution(Multivariate, PointMass, m=[7.0, 9.0]))
end

end # module
