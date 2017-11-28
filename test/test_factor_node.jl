module FactorNodeTest

using Base.Test
import ForneyLab: FactorGraph, FactorNode, Clamp, Variable, Interface, PointMass, GaussianMixture

@testset "FactorNode" begin
    g = FactorGraph()

    # Collect all node types
    node_types = DataType[]
    stack = [FactorNode]
    while !isempty(stack)
        for node_type in subtypes(pop!(stack))
            if isleaftype(node_type)
                (node_type == PointMass) && continue # skip PointMass
                push!(node_types, node_type)
            else
                push!(stack, node_type)
            end
        end
    end

    for node_type in node_types
        # Instantiate a test node
        if node_type == Clamp
            test_node = Clamp(Variable(), 0.0)
        elseif node_type == GaussianMixture # Required for Vararg argument
            test_node = GaussianMixture(Variable(), Variable(), Variable(), Variable(), Variable(), Variable())
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
        if isa(node_type, Clamp)
            @test_throws Exception Clamp(Variable(), 0.0; id=test_node.id)
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
    @test isless(Clamp(Variable(), 0.0; id=:a), Clamp(Variable(), 0.0; id=:b))
end

end #module
