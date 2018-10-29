module VariableTest

using Test
using ForneyLab

@testset "Variable" begin
    g = FactorGraph()

    # Variable should construct and be assigned to graph
    var = Variable(id=:my_var)
    @test isa(var, Variable)
    @test length(g.variables) == 1
    @test ===(g.variables[:my_var], var)
end

@testset "associate!" begin
    FactorGraph()
    var = Variable()

    # Variable should be associated with a new edge
    node1 = Clamp(var, 0.0)
    iface1 = node1.interfaces[1]
    edge1 = iface1.edge
    @test length(var.edges) == 1 # Check Variable
    @test ===(var.edges[1], edge1)
    @test ===(edge1.variable, var) # Check Edge
    @test ===(edge1.a, iface1)
    @test edge1.b == nothing
    @test ===(iface1.edge, edge1) # Check Interface
    @test iface1.partner == nothing

    # Variable should be associated with an existing edge
    node2 = Clamp(var, 0.0)
    iface2 = node2.interfaces[1]
    @test length(var.edges) == 1 # Check Variable
    @test ===(var.edges[1], edge1)
    @test ===(edge1.variable, var) # Check edge
    @test ===(edge1.a, iface1)
    @test ===(edge1.b, iface2)
    @test ===(iface1.edge, edge1) # Check Interfaces
    @test ===(iface2.edge, edge1)
    @test ===(iface1.partner, iface2)
    @test ===(iface2.partner, iface1)

    # Variable should become equality constrained
    node3 = Clamp(var, 0.0)
    eq = node1.interfaces[1].partner.node
    @test isa(eq, Equality) # Check new Equality
    eq_iface1 = eq.interfaces[1]
    eq_iface2 = eq.interfaces[2]
    eq_iface3 = eq.interfaces[3]
    iface3 = node3.interfaces[1]
    edge2 = iface2.edge
    edge3 = iface3.edge
    @test length(var.edges) == 3 # Check Variable
    @test ===(var.edges[1], edge1)
    @test ===(var.edges[2], edge2)
    @test ===(var.edges[3], edge3)
    @test ===(edge1.variable, var) # Check Edges
    @test ===(edge2.variable, var)
    @test ===(edge3.variable, var)
    @test ===(edge1.a, iface1)
    @test ===(edge1.b, eq_iface1)
    @test ===(edge2.a, iface2)
    @test ===(edge2.b, eq_iface2)
    @test ===(edge3.a, iface3)
    @test ===(edge3.b, eq_iface3)
    @test ===(iface1.edge, edge1) # Check Interfaces
    @test ===(eq_iface1.edge, edge1)
    @test ===(iface1.partner, eq_iface1)
    @test ===(eq_iface1.partner, iface1)
    @test ===(iface2.edge, edge2)
    @test ===(eq_iface2.edge, edge2)
    @test ===(iface2.partner, eq_iface2)
    @test ===(eq_iface2.partner, iface2)
    @test ===(iface3.edge, edge3)
    @test ===(eq_iface3.edge, edge3)
    @test ===(iface3.partner, eq_iface3)
    @test ===(eq_iface3.partner, iface3)
end

end #module