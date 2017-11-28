module EqualityTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPEqualityGaussian, SPEqualityGammaWishart, SPEqualityBernoulli, SPEqualityBeta, SPEqualityCategorical, SPEqualityDirichlet, SPEqualityPointMass


#-------------
# Update rules
#-------------

@testset "SPEqualityGaussian" begin
    @test SPEqualityGaussian <: SumProductRule{Equality}
    @test outboundType(SPEqualityGaussian) == Message{Gaussian}
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Message{Gaussian}, Void]) 
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Void, Message{Gaussian}]) 
    @test isApplicable(SPEqualityGaussian, [Void, Message{Gaussian}, Message{Gaussian}]) 

    @test ruleSPEqualityGaussian(Message(Univariate, Gaussian, xi=1.0, w=2.0), Message(Univariate, Gaussian, xi=3.0, w=4.0), nothing) == Message(Univariate, Gaussian, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(Message(Univariate, Gaussian, xi=1.0, w=2.0), nothing, Message(Univariate, Gaussian, xi=3.0, w=4.0)) == Message(Univariate, Gaussian, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(nothing, Message(Univariate, Gaussian, xi=1.0, w=2.0), Message(Univariate, Gaussian, xi=3.0, w=4.0)) == Message(Univariate, Gaussian, xi=4.0, w=6.0)

    @test ruleSPEqualityGaussian(Message(Multivariate, Gaussian, xi=[1.0], w=[2.0].'), Message(Multivariate, Gaussian, xi=[3.0], w=[4.0].'), nothing) == Message(Multivariate, Gaussian, xi=[4.0], w=[6.0].')
    @test ruleSPEqualityGaussian(Message(Multivariate, Gaussian, xi=[1.0], w=[2.0].'), nothing, Message(Multivariate, Gaussian, xi=[3.0], w=[4.0].')) == Message(Multivariate, Gaussian, xi=[4.0], w=[6.0].')
    @test ruleSPEqualityGaussian(nothing, Message(Multivariate, Gaussian, xi=[1.0], w=[2.0].'), Message(Multivariate, Gaussian, xi=[3.0], w=[4.0].')) == Message(Multivariate, Gaussian, xi=[4.0], w=[6.0].')
end

@testset "SPEqualityGammaWishart" begin
    @test SPEqualityGammaWishart <: SumProductRule{Equality}
    @test outboundType(SPEqualityGammaWishart) == Message{Union{Gamma, Wishart}}
    @test isApplicable(SPEqualityGammaWishart, [Message{Union{Gamma, Wishart}}, Message{Union{Gamma, Wishart}}, Void]) 
    @test isApplicable(SPEqualityGammaWishart, [Message{Union{Gamma, Wishart}}, Void, Message{Union{Gamma, Wishart}}]) 
    @test isApplicable(SPEqualityGammaWishart, [Void, Message{Union{Gamma, Wishart}}, Message{Union{Gamma, Wishart}}]) 
    @test isApplicable(SPEqualityGammaWishart, [Message{Gamma}, Message{Gamma}, Void]) 
    @test isApplicable(SPEqualityGammaWishart, [Message{Wishart}, Message{Wishart}, Void]) 

    @test ruleSPEqualityGammaWishart(Message(Univariate, Gamma, a=1.0, b=2.0), Message(Univariate, Gamma, a=3.0, b=4.0), nothing) == Message(Univariate, Gamma, a=3.0, b=6.0)
    @test ruleSPEqualityGammaWishart(Message(Univariate, Gamma, a=1.0, b=2.0), nothing, Message(Univariate, Gamma, a=3.0, b=4.0)) == Message(Univariate, Gamma, a=3.0, b=6.0)
    @test ruleSPEqualityGammaWishart(nothing, Message(Univariate, Gamma, a=1.0, b=2.0), Message(Univariate, Gamma, a=3.0, b=4.0)) == Message(Univariate, Gamma, a=3.0, b=6.0)

    @test ruleSPEqualityGammaWishart(Message(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), Message(MatrixVariate, Wishart, nu=6.0, v=[0.125].'), nothing) == Message(MatrixVariate, Wishart, nu=6.0, v=[0.08333333333333336].')
    @test ruleSPEqualityGammaWishart(Message(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), nothing, Message(MatrixVariate, Wishart, nu=6.0, v=[0.125].')) == Message(MatrixVariate, Wishart, nu=6.0, v=[0.08333333333333336].')
    @test ruleSPEqualityGammaWishart(nothing, Message(MatrixVariate, Wishart, nu=2.0, v=[0.25].'), Message(MatrixVariate, Wishart, nu=6.0, v=[0.125].')) == Message(MatrixVariate, Wishart, nu=6.0, v=[0.08333333333333336].')
end

@testset "SPEqualityBernoulli" begin
    @test SPEqualityBernoulli <: SumProductRule{Equality}
    @test outboundType(SPEqualityBernoulli) == Message{Bernoulli}
    @test isApplicable(SPEqualityBernoulli, [Message{Bernoulli}, Message{Bernoulli}, Void]) 
    @test isApplicable(SPEqualityBernoulli, [Message{Bernoulli}, Void, Message{Bernoulli}]) 
    @test isApplicable(SPEqualityBernoulli, [Void, Message{Bernoulli}, Message{Bernoulli}]) 

    @test ruleSPEqualityBernoulli(Message(Univariate, Bernoulli, p=0.2), Message(Univariate, Bernoulli, p=0.8), nothing) == Message(Univariate, Bernoulli, p=0.5000000000000001)
    @test ruleSPEqualityBernoulli(Message(Univariate, Bernoulli, p=0.2), nothing, Message(Univariate, Bernoulli, p=0.8)) == Message(Univariate, Bernoulli, p=0.5000000000000001)
    @test ruleSPEqualityBernoulli(nothing, Message(Univariate, Bernoulli, p=0.2), Message(Univariate, Bernoulli, p=0.8)) == Message(Univariate, Bernoulli, p=0.5000000000000001)
end

@testset "SPEqualityBeta" begin
    @test SPEqualityBeta <: SumProductRule{Equality}
    @test outboundType(SPEqualityBeta) == Message{Beta}
    @test isApplicable(SPEqualityBeta, [Message{Beta}, Message{Beta}, Void]) 
    @test isApplicable(SPEqualityBeta, [Message{Beta}, Void, Message{Beta}]) 
    @test isApplicable(SPEqualityBeta, [Void, Message{Beta}, Message{Beta}]) 

    @test ruleSPEqualityBeta(Message(Univariate, Beta, a=1.0, b=2.0), Message(Univariate, Beta, a=3.0, b=4.0), nothing) == Message(Univariate, Beta, a=3.0, b=5.0)
    @test ruleSPEqualityBeta(Message(Univariate, Beta, a=1.0, b=2.0), nothing, Message(Univariate, Beta, a=3.0, b=4.0)) == Message(Univariate, Beta, a=3.0, b=5.0)
    @test ruleSPEqualityBeta(nothing, Message(Univariate, Beta, a=1.0, b=2.0), Message(Univariate, Beta, a=3.0, b=4.0)) == Message(Univariate, Beta, a=3.0, b=5.0)
end

@testset "SPEqualityCategorical" begin
    @test SPEqualityCategorical <: SumProductRule{Equality}
    @test outboundType(SPEqualityCategorical) == Message{Categorical}
    @test isApplicable(SPEqualityCategorical, [Message{Categorical}, Message{Categorical}, Void]) 
    @test isApplicable(SPEqualityCategorical, [Message{Categorical}, Void, Message{Categorical}]) 
    @test isApplicable(SPEqualityCategorical, [Void, Message{Categorical}, Message{Categorical}]) 

    @test ruleSPEqualityCategorical(Message(Univariate, Categorical, p=[0.3, 0.7]), Message(Univariate, Categorical, p=[0.7, 0.3]), nothing) == Message(Univariate, Categorical, p=[0.5, 0.5])
    @test ruleSPEqualityCategorical(Message(Univariate, Categorical, p=[0.3, 0.7]), nothing, Message(Univariate, Categorical, p=[0.7, 0.3])) == Message(Univariate, Categorical, p=[0.5, 0.5])
    @test ruleSPEqualityCategorical(nothing, Message(Univariate, Categorical, p=[0.3, 0.7]), Message(Univariate, Categorical, p=[0.7, 0.3])) == Message(Univariate, Categorical, p=[0.5, 0.5])
end

@testset "SPEqualityDirichlet" begin
    @test SPEqualityDirichlet <: SumProductRule{Equality}
    @test outboundType(SPEqualityDirichlet) == Message{Dirichlet}
    @test isApplicable(SPEqualityDirichlet, [Message{Dirichlet}, Message{Dirichlet}, Void]) 
    @test isApplicable(SPEqualityDirichlet, [Message{Dirichlet}, Void, Message{Dirichlet}]) 
    @test isApplicable(SPEqualityDirichlet, [Void, Message{Dirichlet}, Message{Dirichlet}]) 

    @test ruleSPEqualityDirichlet(Message(Multivariate, Dirichlet, a=[1.0, 2.0]), Message(Multivariate, Dirichlet, a=[3.0, 4.0]), nothing) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
    @test ruleSPEqualityDirichlet(Message(Multivariate, Dirichlet, a=[1.0, 2.0]), nothing, Message(Multivariate, Dirichlet, a=[3.0, 4.0])) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
    @test ruleSPEqualityDirichlet(nothing, Message(Multivariate, Dirichlet, a=[1.0, 2.0]), Message(Multivariate, Dirichlet, a=[3.0, 4.0])) == Message(Multivariate, Dirichlet, a=[3.0, 5.0])
end

@testset "SPEqualityPointMass" begin
    @test SPEqualityPointMass <: SumProductRule{Equality}
    @test outboundType(SPEqualityPointMass) == Message{PointMass}
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{Gaussian}, Void]) 
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Void, Message{Gaussian}]) 
    @test isApplicable(SPEqualityPointMass, [Void, Message{PointMass}, Message{Gaussian}]) 
    @test isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{Union{Gamma, Wishart}}, Void]) 
    @test !isApplicable(SPEqualityPointMass, [Message{PointMass}, Message{PointMass}, Void]) 
    @test !isApplicable(SPEqualityPointMass, [Message{Gaussian}, Message{Gaussian}, Void]) 

    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), Message(Univariate, Gaussian, m=0.0, v=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(Message(Univariate, Gaussian, m=0.0, v=1.0), Message(Univariate, PointMass, m=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), nothing, Message(Univariate, Gaussian, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=1.0)
    @test ruleSPEqualityPointMass(nothing, Message(Univariate, PointMass, m=1.0), Message(Univariate, Gaussian, m=0.0, v=1.0)) == Message(Univariate, PointMass, m=1.0)

    @test ruleSPEqualityPointMass(Message(Multivariate, PointMass, m=[1.0, 2.0]), Message(Multivariate, Gaussian, m=zeros(2), v=diageye(2)), nothing) == Message(Multivariate, PointMass, m=[1.0, 2.0])

    @test ruleSPEqualityPointMass(Message(MatrixVariate, PointMass, m=diageye(2)), Message(MatrixVariate, Wishart, nu=3.0, v=diageye(2)), nothing) == Message(MatrixVariate, PointMass, m=diageye(2))

    @test ruleSPEqualityPointMass(Message(Univariate, PointMass, m=1.0), Message(Univariate, Gamma, a=1.0, b=1.0), nothing) == Message(Univariate, PointMass, m=1.0)
    @test_throws Exception ruleSPEqualityPointMass(Message(Univariate, PointMass, m=-1.0), Message(Univariate, Gamma, a=1.0, b=1.0), nothing)
end

end #module