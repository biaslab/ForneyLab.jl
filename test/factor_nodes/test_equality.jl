module EqualityTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPEqualityGaussian, SPEqualityGammaWishart, SPEqualityBernoulli, SPEqualityBeta, SPEqualityCategorical, SPEqualityDirichlet, SPEqualityPointMass, SPEqualityFn, SPEqualityFnG, SPEqualityFnFactor, SPEqualityGFactor, SPEqualityFactor


#-------------
# Update rules
#-------------

@testset "SPEqualityGaussian" begin
    @test SPEqualityGaussian <: SumProductRule{Equality}
    @test outboundType(SPEqualityGaussian) == Message{GaussianWeightedMeanPrecision}
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Message{Gaussian}, Nothing])
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Nothing, Message{Gaussian}])
    @test isApplicable(SPEqualityGaussian, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleSPEqualityGaussian(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0), Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=4.0), nothing) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0), nothing, Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=4.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(nothing, Message(Univariate, GaussianWeightedMeanPrecision, xi=1.0, w=2.0), Message(Univariate, GaussianWeightedMeanPrecision, xi=3.0, w=4.0)) == Message(Univariate, GaussianWeightedMeanPrecision, xi=4.0, w=6.0)

    @test ruleSPEqualityGaussian(Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0)), Message(Multivariate, GaussianWeightedMeanPrecision, xi=[3.0], w=mat(4.0)), nothing) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[4.0], w=mat(6.0))
    @test ruleSPEqualityGaussian(Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0)), nothing, Message(Multivariate, GaussianWeightedMeanPrecision, xi=[3.0], w=mat(4.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[4.0], w=mat(6.0))
    @test ruleSPEqualityGaussian(nothing, Message(Multivariate, GaussianWeightedMeanPrecision, xi=[1.0], w=mat(2.0)), Message(Multivariate, GaussianWeightedMeanPrecision, xi=[3.0], w=mat(4.0))) == Message(Multivariate, GaussianWeightedMeanPrecision, xi=[4.0], w=mat(6.0))
end

@testset "SPEqualityGammaWishart" begin
    @test SPEqualityGammaWishart <: SumProductRule{Equality}
    @test outboundType(SPEqualityGammaWishart) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SPEqualityGammaWishart, [Message{Union{Gamma, Wishart}}, Message{Union{Gamma, Wishart}}, Nothing])
    @test isApplicable(SPEqualityGammaWishart, [Message{Union{Gamma, Wishart}}, Nothing, Message{Union{Gamma, Wishart}}])
    @test isApplicable(SPEqualityGammaWishart, [Nothing, Message{Union{Gamma, Wishart}}, Message{Union{Gamma, Wishart}}])
    @test isApplicable(SPEqualityGammaWishart, [Message{Gamma}, Message{Gamma}, Nothing])
    @test isApplicable(SPEqualityGammaWishart, [Message{Wishart}, Message{Wishart}, Nothing])

    @test ruleSPEqualityGammaWishart(Message(Univariate, Gamma, a=1.0, b=2.0), Message(Univariate, Gamma, a=3.0, b=4.0), nothing) == Message(Univariate, Gamma, a=3.0, b=6.0)
    @test ruleSPEqualityGammaWishart(Message(Univariate, Gamma, a=1.0, b=2.0), nothing, Message(Univariate, Gamma, a=3.0, b=4.0)) == Message(Univariate, Gamma, a=3.0, b=6.0)
    @test ruleSPEqualityGammaWishart(nothing, Message(Univariate, Gamma, a=1.0, b=2.0), Message(Univariate, Gamma, a=3.0, b=4.0)) == Message(Univariate, Gamma, a=3.0, b=6.0)

    @test ruleSPEqualityGammaWishart(Message(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.125)), nothing) == Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.08333333333333336))
    @test ruleSPEqualityGammaWishart(Message(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), nothing, Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.125))) == Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.08333333333333336))
    @test ruleSPEqualityGammaWishart(nothing, Message(MatrixVariate, Wishart, nu=2.0, v=mat(0.25)), Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.125))) == Message(MatrixVariate, Wishart, nu=6.0, v=mat(0.08333333333333336))
end

@testset "SPEqualityBernoulli" begin
    @test SPEqualityBernoulli <: SumProductRule{Equality}
    @test outboundType(SPEqualityBernoulli) == Message{Bernoulli}
    @test isApplicable(SPEqualityBernoulli, [Message{Bernoulli}, Message{Bernoulli}, Nothing])
    @test isApplicable(SPEqualityBernoulli, [Message{Bernoulli}, Nothing, Message{Bernoulli}])
    @test isApplicable(SPEqualityBernoulli, [Nothing, Message{Bernoulli}, Message{Bernoulli}])

    @test ruleSPEqualityBernoulli(Message(Univariate, Bernoulli, p=0.2), Message(Univariate, Bernoulli, p=0.8), nothing) == Message(Univariate, Bernoulli, p=0.5000000000000001)
    @test ruleSPEqualityBernoulli(Message(Univariate, Bernoulli, p=0.2), nothing, Message(Univariate, Bernoulli, p=0.8)) == Message(Univariate, Bernoulli, p=0.5000000000000001)
    @test ruleSPEqualityBernoulli(nothing, Message(Univariate, Bernoulli, p=0.2), Message(Univariate, Bernoulli, p=0.8)) == Message(Univariate, Bernoulli, p=0.5000000000000001)
end

@testset "SPEqualityBeta" begin
    @test SPEqualityBeta <: SumProductRule{Equality}
    @test outboundType(SPEqualityBeta) == Message{Beta}
    @test isApplicable(SPEqualityBeta, [Message{Beta}, Message{Beta}, Nothing])
    @test isApplicable(SPEqualityBeta, [Message{Beta}, Nothing, Message{Beta}])
    @test isApplicable(SPEqualityBeta, [Nothing, Message{Beta}, Message{Beta}])

    @test ruleSPEqualityBeta(Message(Univariate, Beta, a=1.0, b=2.0), Message(Univariate, Beta, a=3.0, b=4.0), nothing) == Message(Univariate, Beta, a=3.0, b=5.0)
    @test ruleSPEqualityBeta(Message(Univariate, Beta, a=1.0, b=2.0), nothing, Message(Univariate, Beta, a=3.0, b=4.0)) == Message(Univariate, Beta, a=3.0, b=5.0)
    @test ruleSPEqualityBeta(nothing, Message(Univariate, Beta, a=1.0, b=2.0), Message(Univariate, Beta, a=3.0, b=4.0)) == Message(Univariate, Beta, a=3.0, b=5.0)
end

@testset "SPEqualityCategorical" begin
    @test SPEqualityCategorical <: SumProductRule{Equality}
    @test outboundType(SPEqualityCategorical) == Message{Categorical}
    @test isApplicable(SPEqualityCategorical, [Message{Categorical}, Message{Categorical}, Nothing])
    @test isApplicable(SPEqualityCategorical, [Message{Categorical}, Nothing, Message{Categorical}])
    @test isApplicable(SPEqualityCategorical, [Nothing, Message{Categorical}, Message{Categorical}])

    @test ruleSPEqualityCategorical(Message(Univariate, Categorical, p=[0.3, 0.7]), Message(Univariate, Categorical, p=[0.7, 0.3]), nothing) == Message(Univariate, Categorical, p=[0.5, 0.5])
    @test ruleSPEqualityCategorical(Message(Univariate, Categorical, p=[0.3, 0.7]), nothing, Message(Univariate, Categorical, p=[0.7, 0.3])) == Message(Univariate, Categorical, p=[0.5, 0.5])
    @test ruleSPEqualityCategorical(nothing, Message(Univariate, Categorical, p=[0.3, 0.7]), Message(Univariate, Categorical, p=[0.7, 0.3])) == Message(Univariate, Categorical, p=[0.5, 0.5])
end

@testset "SPEqualityDirichlet" begin
    @test SPEqualityDirichlet <: SumProductRule{Equality}
    @test outboundType(SPEqualityDirichlet) == Message{Dirichlet}
    @test isApplicable(SPEqualityDirichlet, [Message{Dirichlet}, Message{Dirichlet}, Nothing])
    @test isApplicable(SPEqualityDirichlet, [Message{Dirichlet}, Nothing, Message{Dirichlet}])
    @test isApplicable(SPEqualityDirichlet, [Nothing, Message{Dirichlet}, Message{Dirichlet}])

    @test ruleSPEqualityDirichlet(Message(Multivariate, Dirichlet, a=[1.0, 2.0]), Message(Multivariate, Dirichlet, a=[3.0, 4.0]), nothing) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
    @test ruleSPEqualityDirichlet(Message(Multivariate, Dirichlet, a=[1.0, 2.0]), nothing, Message(Multivariate, Dirichlet, a=[3.0, 4.0])) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
    @test ruleSPEqualityDirichlet(nothing, Message(Multivariate, Dirichlet, a=[1.0, 2.0]), Message(Multivariate, Dirichlet, a=[3.0, 4.0])) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
end

@testset "SPEqualityPointMass" begin
    @test SPEqualityPointMass <: SumProductRule{Equality}
    @test outboundType(SPEqualityPointMass) == Message{PointMass}
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{Gaussian}, Nothing])
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Nothing, Message{Gaussian}])
    @test isApplicable(SPEqualityPointMass, [Nothing, Message{PointMass}, Message{Gaussian}])
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{Union{Gamma, Wishart}}, Nothing])
    @test !isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{PointMass}, Nothing])
    @test !isApplicable(SPEqualityPointMass, [Message{Gaussian}, Message{Gaussian}, Nothing])

    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=1.0)

    @test ruleSPEqualityPointMass(Message(Multivariate, PointMass, m=[1.0, 2.0]), Message(Multivariate, GaussianMeanVariance, m=zeros(2), v=diageye(2)), nothing) == Message(Multivariate, PointMass, m=[1.0, 2.0])

    @test ruleSPEqualityPointMass(Message(MatrixVariate, PointMass, m=diageye(2)), Message(MatrixVariate, Wishart, nu=3.0, v=diageye(2)), nothing) == Message(MatrixVariate, PointMass, m=diageye(2))

    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), Message(Univariate, Gamma, a=1.0, b=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test_throws Exception ruleSPEqualityPointMass(Message(Univariate, PointMass, m=-1.0), Message(Univariate, Gamma, a=1.0, b=1.0), nothing)
end

@testset "SPEqualityFn" begin
    @test SPEqualityFn <: SumProductRule{Equality}
    @test outboundType(SPEqualityFn) == Message{Function}
    @test isApplicable(SPEqualityFn, [Message{Function}, Message{Function}, Nothing])
    @test !isApplicable(SPEqualityFn, [Message{Beta}, Nothing, Message{Gamma}])
    @test !isApplicable(SPEqualityFn, [Nothing, Message{Beta}, Message{Gaussian}])
end

@testset "SPEqualityFnG" begin
    @test SPEqualityFnG <: SumProductRule{Equality}
    @test outboundType(SPEqualityFnG) == Message{GaussianMeanVariance}
    @test isApplicable(SPEqualityFnG, [Nothing, Message{Function}, Message{Gaussian}])
    @test !isApplicable(SPEqualityFnG, [Message{Function}, Message{Function}, Nothing])
    @test !isApplicable(SPEqualityFnG, [Message{Beta}, Nothing, Message{Gamma}])
    @test !isApplicable(SPEqualityFnG, [Nothing, Message{Beta}, Message{Gaussian}])
end

@testset "SPEqualityFnFactor" begin
    @test SPEqualityFnFactor <: SumProductRule{Equality}
    @test outboundType(SPEqualityFnFactor) == Message{SampleList}
    @test isApplicable(SPEqualityFnFactor, [Nothing, Message{Function}, Message{Poisson}])
    @test !isApplicable(SPEqualityFnFactor, [Nothing, Message{Function}, Message{Gaussian}])
    @test !isApplicable(SPEqualityFnFactor, [Message{Function}, Message{Function}, Nothing])
    @test !isApplicable(SPEqualityFnFactor, [Message{Beta}, Nothing, Message{Gamma}])
    @test !isApplicable(SPEqualityFnFactor, [Nothing, Message{Beta}, Message{Gaussian}])
end

@testset "SPEqualityGFactor" begin
    @test SPEqualityGFactor <: SumProductRule{Equality}
    @test outboundType(SPEqualityGFactor) == Message{GaussianMeanVariance}
    @test !isApplicable(SPEqualityGFactor, [Nothing, Message{Function}, Message{Poisson}])
    @test !isApplicable(SPEqualityGFactor, [Nothing, Message{Function}, Message{Gaussian}])
    @test !isApplicable(SPEqualityGFactor, [Message{Function}, Message{Function}, Nothing])
    @test !isApplicable(SPEqualityGFactor, [Message{Beta}, Nothing, Message{Gamma}])
    @test isApplicable(SPEqualityGFactor, [Nothing, Message{Beta}, Message{Gaussian}])
end

@testset "SPEqualityFactor" begin
    @test SPEqualityFactor <: SumProductRule{Equality}
    @test outboundType(SPEqualityFactor) == Message{SampleList}
    @test !isApplicable(SPEqualityFactor, [Nothing, Message{Function}, Message{Poisson}])
    @test !isApplicable(SPEqualityFactor, [Nothing, Message{Function}, Message{Gaussian}])
    @test !isApplicable(SPEqualityFactor, [Message{Function}, Message{Function}, Nothing])
    @test isApplicable(SPEqualityFactor, [Message{Beta}, Nothing, Message{Gamma}])
    @test !isApplicable(SPEqualityFactor, [Message{Beta}, Nothing, Message{Beta}])
    @test !isApplicable(SPEqualityFactor, [Nothing, Message{Beta}, Message{Gaussian}])
end

end #module
