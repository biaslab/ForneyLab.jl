module InterfaceTest

using Test
import ForneyLab: Interface, FactorNode

# Integration helper
mutable struct MockNode <: FactorNode end

@testset "Interface" begin
    # Interface should construct
    node = MockNode()
    iface = Interface(node)
    @test ===(iface.node, node)
    @test iface.edge == nothing
    @test iface.partner == nothing
end

end # module