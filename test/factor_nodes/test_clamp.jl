module ClampTest

using Base.Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, SPClamp

@testset "Clamp" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)

    @test isa(nd, Clamp)
    @test nd.value == 1.0
end

@testset "constant" begin
    g = FactorGraph()
    var = constant(1.0, id=:my_constant)
    nd = g.nodes[:my_constant]

    @test isa(var, Variable)
    @test isa(nd, Clamp)
    @test nd.value == 1.0
end

@testset "placeholder" begin
    g = FactorGraph()

    # Standard placeholder
    var = Variable()
    placeholder(var, :y)
    nd = g.nodes[:placeholder_y]

    @test isa(nd, Clamp)
    @test g.placeholders[nd] == (:y, 0)

    # Indexed placeholder
    var_i = Variable()
    placeholder(var_i, :y, index=1)
    nd_i = g.nodes[:placeholder_y_1]

    @test isa(nd_i, Clamp)
    @test g.placeholders[nd_i] == (:y, 1)
end


#-------------
# Update rules
#-------------

@testset "SPClamp" begin
    @test SPClamp <: SumProductRule{Clamp}
    @test outboundType(SPClamp) == Message{PointMass}
    @test isApplicable(SPClamp, DataType[]) 
end

end #module