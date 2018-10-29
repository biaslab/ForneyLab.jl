module DependencyGraphTest

using Test
import ForneyLab: LinkedList, DependencyGraph, addVertex!, addEdge!, children, neighbors

mutable struct MockVertex
    id::Symbol
end

@testset "LinkedList" begin
    ll = LinkedList{typeof("test")}()
    @test isempty(ll) == true
    @test_throws Exception push!(ll, 3) # Incorrect element type
    push!(ll, "test1")
    push!(ll, "test2")
    @test isempty(ll) == false
    @test collect(ll) == ["test1"; "test2"]
end

@testset "DependencyGraph" begin
    @testset "DependencyGraph()" begin
        # should construct a dependency graph
        dg = DependencyGraph{MockVertex}()
        @test isa(dg, DependencyGraph)
        @test dg.vertices == MockVertex[]
        @test dg.neighbors == LinkedList[]

        # Construct:
        # (1)-->(2)

        # addVertex! should add a vertex
        v1 = MockVertex(:v1)
        v2 = MockVertex(:v2)
        addVertex!(dg, v1)
        addVertex!(dg, v2)
        @test dg.vertices == MockVertex[v1, v2]

        # addEdge! should construct an edge
        addEdge!(dg, v1, v2)
        @test !(v1 in neighbors(v2, dg))
        @test (v2 in neighbors(v1, dg)) # Neighbours is not a commutative relation
        @test !(v1 in neighbors(v1, dg)) # A vertex is not a neighbor of itself
        @test !(v2 in neighbors(v2, dg))

        # Extend:
        # (1)-->(2)-->(3)
        v3 = MockVertex(:v3)
        addVertex!(dg, v3)
        addEdge!(dg, v2, v3)
        @test !(v2 in neighbors(v3, dg))
        @test (v3 in neighbors(v2, dg))
        @test !(v1 in neighbors(v3, dg)) # Unconnected vertices are not neighbours
        @test !(v3 in neighbors(v1, dg))
    end

    @testset "children" begin
        # should find all children of a vertex and return a topologically sorted array

        # Construct:
        #
        #  (1)-->(4)<--(5)
        #   |     |     |
        #   v     v     v
        #  (2)   (6)   (7)
        #   |
        #   v
        #  (3)

        dg = DependencyGraph{MockVertex}()
        vertices = MockVertex[]
        for i = 1:7
            v = MockVertex(Symbol("v$i"))
            push!(vertices, v)
            addVertex!(dg, v)
        end
        addEdge!(dg, vertices[1], vertices[2])
        addEdge!(dg, vertices[2], vertices[3])
        addEdge!(dg, vertices[1], vertices[4])
        addEdge!(dg, vertices[4], vertices[6])
        addEdge!(dg, vertices[5], vertices[4])
        addEdge!(dg, vertices[5], vertices[7])

        @test children(vertices[1], dg) == MockVertex[vertices[3], vertices[2], vertices[6], vertices[4], vertices[1]]
        @test children(vertices[1], dg, breaker_sites=Set{MockVertex}([vertices[4]])) == MockVertex[vertices[3], vertices[2], vertices[1]]
        @test children(vertices[1], dg, restrict_to=Set{MockVertex}(vertices[1:4])) == MockVertex[vertices[3], vertices[2], vertices[4], vertices[1]]
    end

    @testset "Loopy dependengy graph" begin
        # should find all children of a vertex and return a topologically sorted array in a cyclic graph

        # Construct:
        #
        #  (1)-->(2)
        #   ^    /
        #    \  v
        #     (3)

        dg = DependencyGraph{MockVertex}()
        vertices = MockVertex[]
        for i = 1:3
            v = MockVertex(Symbol("v$i"))
            push!(vertices, v)
            addVertex!(dg, v)
        end
        addEdge!(dg, vertices[1], vertices[2])
        addEdge!(dg, vertices[2], vertices[3])
        addEdge!(dg, vertices[3], vertices[1])

        @test (vertices[1] in neighbors(vertices[3], dg))
        @test !(vertices[2] in neighbors(vertices[3], dg))
        @test (vertices[2] in neighbors(vertices[1], dg))
        @test !(vertices[3] in neighbors(vertices[1], dg))
        @test (vertices[3] in neighbors(vertices[2], dg))
        @test !(vertices[1] in neighbors(vertices[2], dg))

        @test_throws Exception children(vertices[1], dg) # children function errors when cycles are not allowed
        @test children(vertices[1], dg, allow_cycles=true) == MockVertex[vertices[3], vertices[2], vertices[1]]
    end
end

end # module