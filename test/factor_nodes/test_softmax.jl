module SoftmaxTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: VBSoftmaxOut, VBSoftmaxIn1, VBSoftmaxXi, VBSoftmaxA


#-------------
# Update rules
#-------------

@testset "VBSoftmaxOut" begin
    @test VBSoftmaxOut <: NaiveVariationalRule{Softmax}
    @test outboundType(VBSoftmaxOut) == Message{Categorical}
    @test isApplicable(VBSoftmaxOut, [Nothing, Distribution, Distribution, Distribution]) 

    @test ruleVBSoftmaxOut(nothing, Distribution(Multivariate, Gaussian{Moments}, m=[-1.0, 1.0], v=diageye(2)), Distribution(Multivariate, PointMass, m=[1.0, 2.0]), Distribution(Univariate, PointMass, m=1.0)) == Message(Univariate, Categorical, p=[0.11920292202211755, 0.8807970779778824])
end

@testset "VBSoftmaxIn1" begin
    @test VBSoftmaxIn1 <: NaiveVariationalRule{Softmax}
    @test outboundType(VBSoftmaxIn1) == Message{Gaussian{Canonical}}
    @test isApplicable(VBSoftmaxIn1, [Distribution, Nothing, Distribution, Distribution]) 

    @test ruleVBSoftmaxIn1(Distribution(Univariate, Categorical, p=[0.8, 0.2]), nothing, Distribution(Multivariate, PointMass, m=[1.0, 2.0]), Distribution(Univariate, PointMass, m=1.0)) == Message(Multivariate, Gaussian{Canonical}, xi=[0.5097033519265977, -0.039212315743835824], w=[0.2310585786300049 0.0; 0.0 0.19039853898894119])
end

@testset "VBSoftmaxXi" begin
    @test VBSoftmaxXi <: NaiveVariationalRule{Softmax}
    @test outboundType(VBSoftmaxXi) == Message{Function}
    @test isApplicable(VBSoftmaxXi, [Distribution, Distribution, Nothing, Distribution]) 

    @test ruleVBSoftmaxXi(Distribution(Univariate, Categorical, p=[0.8, 0.2]), Distribution(Multivariate, Gaussian{Moments}, m=[-1.0, 1.0], v=diageye(2)), nothing, Distribution(Univariate, PointMass, m=1.0)) == Message(Multivariate, Function, mode=[2.23606797749979, 1.0])
end

@testset "VBSoftmaxA" begin
    @test VBSoftmaxA <: NaiveVariationalRule{Softmax}
    @test outboundType(VBSoftmaxA) == Message{Function}
    @test isApplicable(VBSoftmaxA, [Distribution, Distribution, Distribution, Nothing]) 

    @test ruleVBSoftmaxA(Distribution(Univariate, Bernoulli, p=[0.8, 0.2]), Distribution(Multivariate, Gaussian{Moments}, m=[-1.0, 1.0], v=diageye(2)), Distribution(Multivariate, PointMass, m=[1.0, 2.0]), nothing) == Message(Univariate, Function, mode=-0.09647491510115125)
end

@testset "averageEnergy" begin
    @test averageEnergy(Softmax, Distribution(Univariate, Categorical, p=[0.8, 0.2]), Distribution(Multivariate, Gaussian{Moments}, m=[-1.0, 1.0], v=diageye(2)), Distribution(Multivariate, PointMass, m=[1.0, 2.0]), Distribution(Univariate, PointMass, m=1.0)) == 2.907107586326735
end

end # module