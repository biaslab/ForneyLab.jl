module VariableTest

using Base.Test
import ForneyLab: Variable, Constant, Equality, FactorGraph

@testset "Variable" begin
    # Variable should construct
    @test isa(Variable(), Variable)
end

@testset "associate!" begin
    FactorGraph()
    var = Variable()

    # Variable should be associated with a new edge
    node1 = Constant(var, 0.0)
    iface1 = node1.interfaces[1]
    edge1 = iface1.edge
    @test length(var.edges) == 1 # Check Variable
    @test is(var.edges[1], edge1)
    @test is(edge1.variable, var) # Check Edge
    @test is(edge1.a, iface1)
    @test edge1.b == nothing    
    @test is(iface1.edge, edge1) # Check Interface
    @test iface1.partner == nothing    

    # Variable should be associated with an existing edge 
    node2 = Constant(var, 0.0)
    iface2 = node2.interfaces[1]
    @test length(var.edges) == 1 # Check Variable
    @test is(var.edges[1], edge1)
    @test is(edge1.variable, var) # Check edge
    @test is(edge1.a, iface1)
    @test is(edge1.b, iface2)
    @test is(iface1.edge, edge1) # Check Interfaces
    @test is(iface2.edge, edge1)
    @test is(iface1.partner, iface2)
    @test is(iface2.partner, iface1)

    # Variable should become equality constrained
    node3 = Constant(var, 0.0)
    eq = node1.interfaces[1].partner.node
    @test isa(eq, Equality) # Check new Equality
    eq_iface1 = eq.interfaces[1]
    eq_iface2 = eq.interfaces[2]
    eq_iface3 = eq.interfaces[3]
    iface3 = node3.interfaces[1]
    edge2 = iface2.edge
    edge3 = iface3.edge
    @test length(var.edges) == 3 # Check Variable
    @test is(var.edges[1], edge1)
    @test is(var.edges[2], edge2)
    @test is(var.edges[3], edge3)
    @test is(edge1.variable, var) # Check Edges
    @test is(edge2.variable, var)
    @test is(edge3.variable, var)
    @test is(edge1.a, iface1)
    @test is(edge1.b, eq_iface1)
    @test is(edge2.a, iface2)
    @test is(edge2.b, eq_iface2)
    @test is(edge3.a, iface3)
    @test is(edge3.b, eq_iface3)
    @test is(iface1.edge, edge1) # Check Interfaces
    @test is(eq_iface1.edge, edge1)
    @test is(iface1.partner, eq_iface1)
    @test is(eq_iface1.partner, iface1)
    @test is(iface2.edge, edge2)
    @test is(eq_iface2.edge, edge2)
    @test is(iface2.partner, eq_iface2)
    @test is(eq_iface2.partner, iface2)
    @test is(iface3.edge, edge3)
    @test is(eq_iface3.edge, edge3)
    @test is(iface3.partner, eq_iface3)
    @test is(eq_iface3.partner, iface3)
end

end #module