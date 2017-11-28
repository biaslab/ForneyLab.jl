module EdgeTest

using Base.Test
import ForneyLab: Interface, Edge, Variable, Interface, FactorNode, FactorGraph, currentGraph, addNode!, disconnect!, generateId

# Integration helper
type MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function MockNode(; id=generateId(MockNode))
        self = new(id, Interface[], Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        return self
    end
end

@testset "Edge" begin
    # Edge should couple interfaces
    FactorGraph()

    a = Interface(MockNode())
    b = Interface(MockNode())
    edge = Edge(Variable(), a, b)

    @test is(edge.a, a)
    @test is(edge.b, b)
    @test is(a.edge, edge)
    @test is(b.edge, edge)
    @test is(a.partner, b)
    @test is(b.partner, a)
end

@testset "disconnect!" begin
    FactorGraph()
    a = Interface(MockNode())
    b = Interface(MockNode())
    edge = Edge(Variable(), a, b)

    # disconnect! from 'a' side should decouple interfaces
    disconnect!(edge, a)
    @test is(edge.a, b)
    @test edge.b == nothing
    @test a.edge == nothing
    @test is(b.edge, edge)
    @test a.partner == nothing
    @test b.partner == nothing

    FactorGraph()
    a = Interface(MockNode())
    b = Interface(MockNode())
    edge = Edge(Variable(), a, b)

    # disconnect! from 'b' side should decouple interfaces
    disconnect!(edge, b)
    @test is(edge.a, a)
    @test edge.b == nothing
    @test is(a.edge, edge)
    @test b.edge == nothing
    @test a.partner == nothing
    @test b.partner == nothing
end

end # module