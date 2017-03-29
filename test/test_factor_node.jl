module FactorNodeTest

using Base.Test
import ForneyLab: FactorGraph, SoftFactor, DeltaFactor, Constant, Variable, Interface, PointMass, Terminal

@testset "FactorNode" begin
    g = FactorGraph()

    for node_type in [subtypes(SoftFactor); subtypes(DeltaFactor); Terminal]
        if node_type == PointMass
            continue
        end

        # Instantiate a test node
        if isa(node_type, Constant)
            test_node = Constant(Variable(), 0.0)
        else
            constructor_argument_length = length(first(methods(node_type)).sig.parameters) - 1
            vars = [Variable() for v = 1:constructor_argument_length]
            test_node = node_type(vars...)
        end

        # Node fields should be of correct types
        @test isa(test_node.id, Symbol)
        @test isa(test_node.interfaces, Vector{Interface})
        @test isa(test_node.i, Dict)

        # Node constructor should automatically assign an id and check uniqueness
        @test !isempty(string(test_node.id))
        if isa(node_type, Constant)
            @test_throws Exception Constant(Variable(), 0.0; id=test_node.id)
        else
            @test_throws Exception node_type(vars...; id=test_node.id)
        end

        # Node constructor should assign interfaces to itself
        for iface in test_node.interfaces
            @test is(iface.node, test_node)
        end

        # Node constructor should add node to graph
        @test is(g.nodes[test_node.id], test_node)

    end

    # isless should be defined for nodes
    @test isless(Constant(Variable(), 0.0; id=:a), Constant(Variable(), 0.0; id=:b))
end

end #module
