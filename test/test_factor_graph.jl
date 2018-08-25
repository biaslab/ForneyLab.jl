module FactorGraphTest

using Test
import ForneyLab: FactorGraph, FactorNode, Interface, Edge, Variable, generateId, currentGraph, setCurrentGraph,
generateId, addNode!, hasNode, addVariable!, hasVariable, Clamp

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function MockNode(; id=generateId(MockNode))
        self = new(id, Interface[], Dict{Symbol,Interface}())

        return self
    end
end

@testset "FactorGraph" begin
    # FactorGraph should initiate a new factor graph
    g = FactorGraph()

    @test isa(g.nodes, Dict{Symbol, FactorNode})
    @test isa(g.edges, Vector{Edge})
    @test isa(g.variables, Dict{Symbol, Variable})
    @test isa(g.counters, Dict{String, Int})
    @test isa(g.placeholders, Dict{Clamp, Tuple{Symbol, Int}})

    # currentGraph() should point to the current graph
    @test ===(currentGraph(), g)
    f = FactorGraph()
    @test ===(currentGraph(), f)
    setCurrentGraph(g)
    @test ===(currentGraph(), g)
end

@testset "generateId" begin
    g = FactorGraph()

    # generateId should generate a unique id based on type
    @test generateId(MockNode) == :mocknode_1
    @test generateId(MockNode) == :mocknode_2
    @test g.counters["mocknode"] == 2
end

@testset "addNode!" begin
    g = FactorGraph()
    nd = MockNode(id=:nd)

    # addNode! should register a node with a graph
    @test isempty(g.nodes)
    addNode!(g, nd)
    @test ===(g.nodes[:nd], nd)
    @test hasNode(g, nd)

    # addNode should not add the same node twice
    @test_throws Exception addNode!(g, nd)
end

end