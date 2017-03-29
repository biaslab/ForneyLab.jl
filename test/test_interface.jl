module InterfaceTest

using Base.Test
import ForneyLab: Interface, FactorNode

# Integration helper
type MockNode <: FactorNode end

@testset "Interface" begin
    # Interface should construct
    node = MockNode()
    iface = Interface(node)
    @test is(iface.node, node)
    @test iface.edge == nothing
    @test iface.partner == nothing
end

end # module