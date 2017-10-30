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
    @test outboundType(SPEqualityGaussian) == Message{Univariate{Gaussian}}
    @test isApplicable(SPEqualityGaussian, [Message{Univariate{Gaussian}}, Message{Univariate{Gaussian}}, Void]) 
    @test isApplicable(SPEqualityGaussian, [Message{Univariate{Gaussian}}, Void, Message{Univariate{Gaussian}}]) 
    @test isApplicable(SPEqualityGaussian, [Void, Message{Univariate{Gaussian}}, Message{Univariate{Gaussian}}]) 

    @test ruleSPEqualityGaussian(Message(Univariate(Gaussian, xi=1.0, w=2.0)), Message(Univariate(Gaussian, xi=3.0, w=4.0)), nothing) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
    @test ruleSPEqualityGaussian(Message(Univariate(Gaussian, xi=1.0, w=2.0)), nothing, Message(Univariate(Gaussian, xi=3.0, w=4.0))) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
    @test ruleSPEqualityGaussian(nothing, Message(Univariate(Gaussian, xi=1.0, w=2.0)), Message(Univariate(Gaussian, xi=3.0, w=4.0))) == Message(Univariate(Gaussian, xi=4.0, w=6.0))
end

@testset "SPEqualityGamma" begin
    @test SPEqualityGamma <: SumProductRule{Equality}
    @test outboundType(SPEqualityGamma) == Message{Univariate{Gamma}}
    @test isApplicable(SPEqualityGamma, [Message{Univariate{Gamma}}, Message{Univariate{Gamma}}, Void]) 
    @test isApplicable(SPEqualityGamma, [Message{Univariate{Gamma}}, Void, Message{Univariate{Gamma}}]) 
    @test isApplicable(SPEqualityGamma, [Void, Message{Univariate{Gamma}}, Message{Univariate{Gamma}}]) 

    @test ruleSPEqualityGamma(Message(Univariate(Gamma, a=1.0, b=2.0)), Message(Univariate(Gamma, a=3.0, b=4.0)), nothing) == Message(Univariate(Gamma, a=3.0, b=6.0))
    @test ruleSPEqualityGamma(Message(Univariate(Gamma, a=1.0, b=2.0)), nothing, Message(Univariate(Gamma, a=3.0, b=4.0))) == Message(Univariate(Gamma, a=3.0, b=6.0))
    @test ruleSPEqualityGamma(nothing, Message(Univariate(Gamma, a=1.0, b=2.0)), Message(Univariate(Gamma, a=3.0, b=4.0))) == Message(Univariate(Gamma, a=3.0, b=6.0))
end

end #module