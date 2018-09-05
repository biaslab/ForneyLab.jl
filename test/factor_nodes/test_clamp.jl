module ClampTest

using Test
using ForneyLab
import ForneyLab: outboundType, isApplicable, @ensureVariables, generateId, addNode!, associate!
import ForneyLab: SPClamp

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function MockNode(out; id=generateId(MockNode))
        @ensureVariables(out)
        self = new(id, Array{Interface}(undef, 1), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)

        return self
    end
end

@testset "@ensureVariables" begin
    g = FactorGraph()
    nd_c = MockNode(constant(1.0))
    @test isa(nd_c.i[:out].partner.node, Clamp)

    nd_v = MockNode(1.0)
    @test isa(nd_v.i[:out].partner.node, Clamp)
end

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
    var = Variable(id=:y)
    placeholder(var, :y)
    nd = g.nodes[:placeholder_y]

    @test isa(nd, Clamp)
    @test g.placeholders[nd] == (:y, 0)

    # Placeholder without explicit variable initialization
    g = FactorGraph()

    placeholder(:y)
    nd = g.nodes[:placeholder_y]

    @test length(g.variables) == 1
    @test isa(nd, Clamp)
    @test g.placeholders[nd] == (:y, 0)

    # Indexed placeholder
    g = FactorGraph()

    var_i = Variable(id=:y)
    placeholder(var_i, :y, index=1)
    nd_i = g.nodes[:placeholder_y]

    @test isa(nd_i, Clamp)
    @test g.placeholders[nd_i] == (:y, 1)
end


#-------------
# Update rules
#-------------

@testset "SPClamp" begin
    @test SPClamp{Univariate} <: SumProductRule{Clamp{Univariate}}
    @test outboundType(SPClamp{Univariate}) == Message{PointMass, Univariate}
    @test isApplicable(SPClamp{Univariate}, DataType[])

    @test SPClamp{Multivariate} <: SumProductRule{Clamp{Multivariate}}
    @test outboundType(SPClamp{Multivariate}) == Message{PointMass, Multivariate}
    @test isApplicable(SPClamp{Multivariate}, DataType[])

    @test SPClamp{MatrixVariate} <: SumProductRule{Clamp{MatrixVariate}}
    @test outboundType(SPClamp{MatrixVariate}) == Message{PointMass, MatrixVariate}
    @test isApplicable(SPClamp{MatrixVariate}, DataType[])
end

end #module