module EqualityTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPEqualityGaussian, SPEqualityGamma


#-------------
# Update rules
#-------------

@testset "SPEqualityGaussian" begin
    @test SPEqualityGaussian <: SumProductRule{Equality}
    @test outboundType(SPEqualityGaussian) == Message{Gaussian}
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Message{Gaussian}, Void]) 
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Void, Message{Gaussian}]) 
    @test isApplicable(SPEqualityGaussian, [Void, Message{Gaussian}, Message{Gaussian}]) 

    @test ruleSPEqualityGaussian(Message(Univariate(Gaussian, xi=1.0, w=2.0)), Message(Univariate(Gaussian, xi=3.0, w=4.0)), nothing) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
    @test ruleSPEqualityGaussian(Message(Univariate(Gaussian, xi=1.0, w=2.0)), nothing, Message(Univariate(Gaussian, xi=3.0, w=4.0))) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
    @test ruleSPEqualityGaussian(nothing, Message(Univariate(Gaussian, xi=1.0, w=2.0)), Message(Univariate(Gaussian, xi=3.0, w=4.0))) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
end

@testset "SPEqualityGamma" begin
    @test SPEqualityGamma <: SumProductRule{Equality}
    @test outboundType(SPEqualityGamma) == Message{AbstractGamma}
    @test isApplicable(SPEqualityGamma, [Message{AbstractGamma}, Message{AbstractGamma}, Void]) 
    @test isApplicable(SPEqualityGamma, [Message{AbstractGamma}, Void, Message{AbstractGamma}]) 
    @test isApplicable(SPEqualityGamma, [Void, Message{AbstractGamma}, Message{AbstractGamma}]) 

    @test ruleSPEqualityGamma(Message(Univariate(Gamma, a=1.0, b=2.0)), Message(Univariate(Gamma, a=3.0, b=4.0)), nothing) == Message(Univariate(Gamma, a=3.0, b=6.0))
    @test ruleSPEqualityGamma(Message(Univariate(Gamma, a=1.0, b=2.0)), nothing, Message(Univariate(Gamma, a=3.0, b=4.0))) == Message(Univariate(Gamma, a=3.0, b=6.0))
    @test ruleSPEqualityGamma(nothing, Message(Univariate(Gamma, a=1.0, b=2.0)), Message(Univariate(Gamma, a=3.0, b=4.0))) == Message(Univariate(Gamma, a=3.0, b=6.0))
end

end #module