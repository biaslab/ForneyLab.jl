module DependencyGraphTest

using Base.Test
import ForneyLab: LinkedList

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
    @testset "children" begin
        # should find all children of a vertex and return a topologically sorted array

        # TODO: implement FactorGraph independent tests here
    end
end

end # module