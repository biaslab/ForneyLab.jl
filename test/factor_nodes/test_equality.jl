module EqualityTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPEqualityGaussian


#-------------
# Update rules
#-------------

@testset "SPEqualityGaussian" begin
    @test SPEqualityGaussian <: SumProductRule{Equality}
    @test outboundType(SPEqualityGaussian) == Message{Gaussian}
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Message{Gaussian}, Void]) 
    @test isApplicable(SPEqualityGaussian, [Message{Gaussian}, Void, Message{Gaussian}]) 
    @test isApplicable(SPEqualityGaussian, [Void, Message{Gaussian}, Message{Gaussian}]) 

    @test ruleSPEqualityGaussian(Message(Gaussian, xi=1.0, w=2.0), Message(Gaussian, xi=3.0, w=4.0), nothing) == Message(Gaussian, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(Message(Gaussian, xi=1.0, w=2.0), nothing, Message(Gaussian, xi=3.0, w=4.0)) == Message(Gaussian, xi=4.0, w=6.0)
    @test ruleSPEqualityGaussian(nothing, Message(Gaussian, xi=1.0, w=2.0), Message(Gaussian, xi=3.0, w=4.0)) == Message(Gaussian, xi=4.0, w=6.0)
end

end #module