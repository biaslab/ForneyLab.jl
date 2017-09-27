module BernoulliTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPBernoulliPV


#-------------
# Update rules
#-------------

@testset "SPBernoulliPV" begin
    @test SPBernoulliPV <: SumProductRule{Bernoulli}
    @test outboundType(SPBernoulliPV) == Message{Bernoulli}
    @test isApplicable(SPBernoulliPV, [Message{PointMass}, Void]) 

    @test ruleSPBernoulliPV(Message(PointMass, m=0.2), nothing) == Message(Bernoulli, p=0.2)
end

end # module