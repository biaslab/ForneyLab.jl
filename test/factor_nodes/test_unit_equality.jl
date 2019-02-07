module UnitEqualityTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable
import ForneyLab: SPUnitEqualityOutVGG, SPUnitEqualityIn1GVG, slug


@testset "SPUnitEqualityOutVGG" begin
    @test SPUnitEqualityOutVGG <: SumProductRule{UnitEquality}
    @test outboundType(SPUnitEqualityOutVGG) == Message{GaussianMeanPrecision}
    @test isApplicable(SPUnitEqualityOutVGG, [Nothing, Message{Gaussian}, Message{Gaussian}])

    @test ruleSPUnitEqualityOutVGG(nothing, Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(2.0)), Message(Univariate, GaussianMeanPrecision, m=3.0, w=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0))
    @test ruleSPUnitEqualityOutVGG(nothing, Message(Multivariate, GaussianMeanPrecision, m=[0.0, 0.0], w=[1.0 0.0; 0.0 1.0]), Message(Univariate, GaussianMeanPrecision, m=0.0, w=1.0)) == Message(Multivariate, GaussianMeanPrecision, m=[0.0, 0.0], w=[2.0 0.0; 0.0 1.0])
end

@testset "SPUnitEqualityIn1GVG" begin
    @test SPUnitEqualityIn1GVG <: SumProductRule{UnitEquality}
    @test outboundType(SPUnitEqualityIn1GVG) == Message{GaussianMeanPrecision}
    @test isApplicable(SPUnitEqualityIn1GVG, [Message{Gaussian}, Nothing, Message{Gaussian}])

    @test ruleSPUnitEqualityIn1GVG(Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(2.0)), nothing, Message(Univariate, GaussianMeanPrecision, m=3.0, w=2.0)) == Message(Multivariate, GaussianMeanPrecision, m=[3.0], w=mat(4.0))
end

@testset "slug" begin
    @test slug(UnitEquality) == "unit_eq"
end

end # module
