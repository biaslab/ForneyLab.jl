module VariableTest

using Base.Test
import ForneyLab: Variable, Constant, GaussianMeanVariance, Equality, FactorGraph, @~, currentGraph, constant

@testset "Variable" begin
    g = FactorGraph()

    # Variable should construct and be assigned to graph
    var = Variable(id=:my_var)
    @test isa(var, Variable)
    @test length(g.variables) == 1
    @test is(g.variables[:my_var], var)
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

@testset "@~" begin
    g = FactorGraph()

    # @~ should construct a new variable
    x = constant(0.0)
    y ~ GaussianMeanVariance(x, constant(1.0))
    @test length(g.variables) == 3 # including constants
    @test haskey(g.variables, y.id)
    @test is(g.variables[y.id], y)

    # @~ should reuse existing variable and handle keyword agruments
    y_old = y
    y ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:g_node)
    @test length(g.variables) == 5 # including constants
    @test haskey(g.nodes, :g_node)
    @test is(y, y_old)

    # @~ should handle array element assignments
    g = FactorGraph()
    vars = Vector{Variable}(2)
    vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst1) # new Variable
    @test length(g.variables) == 3 # including constants
    vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst2) # existing Variable
    @test length(g.variables) == 5 # including constants
end

end #module