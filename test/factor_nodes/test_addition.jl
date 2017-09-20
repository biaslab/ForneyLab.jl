module AdditionTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPAdditionGGV

@testset "Addition node construction through + syntax" begin
    g = FactorGraph()

    x ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    y ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    z = x + y

    @test isa(z, Variable)
    @test isa(g.nodes[:addition_1], Addition)
end


#-------------
# Update rules
#-------------

@testset "SPAdditionGGV" begin
    @test SPAdditionGGV <: SumProductRule{Addition}
    @test outboundType(SPAdditionGGV) == Message{Gaussian}
    @test isApplicable(SPAdditionGGV, [Message{Gaussian}, Message{Gaussian}, Void]) 
    @test !isApplicable(SPAdditionGGV, [Message{PointMass}, Message{PointMass}, Void]) 
    @test !isApplicable(SPAdditionGGV, [Void, Message{Gaussian}, Message{Gaussian}]) 
end

# TODO: Add more tests
@testset "ruleSPAdditionGGV" begin
    @test ruleSPAdditionGGV(Message(Gaussian, m=1.0, v=2.0), Message(Gaussian, m=3.0, v=4.0), nothing) == Message(Gaussian, m=4.0, v=6.0)
end

end # module