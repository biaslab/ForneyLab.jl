module VariableTest

using Test
using ForneyLab

@testset "Variable" begin
    g = FactorGraph()

    # Variable call should construct a Variable that is assigned to graph
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

@testset "@RV" begin
    g = FactorGraph()

    # @RV should construct a new variable
    x = constant(0.0)
    @RV y ~ GaussianMeanVariance(x, constant(1.0))
    @test length(g.variables) == 3 # including constants
    @test y.id == :y # automatically assigned id based on variable name in code
    @test haskey(g.variables, y.id)
    @test g.variables[y.id] === y

    # @RV ~ should reuse existing variable and handle keyword agruments
    y_old = y
    @RV y ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:g_node)
    @test length(g.variables) == 5 # including constants
    @test haskey(g.nodes, :g_node)
    @test y === y_old

    # @RV should handle array element assignments and explicit Variable ids
    g = FactorGraph()
    vars = Vector{Variable}(undef, 2)
    i = 1
    @RV [id=:v*i] vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst1) # new Variable
    @test length(g.variables) == 3 # including constants
    @test vars[1].id == :v1
    @RV vars[1] ~ GaussianMeanVariance(constant(0.0), constant(1.0); id=:tst2) # existing Variable
    @test length(g.variables) == 5 # including constants
    @RV vars[2*i] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test vars[2*i].id == :vars_2
    varmatrix = Matrix{Variable}(undef,2,2)
    @RV varmatrix[1,2*i] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test varmatrix[1,2*i].id == :varmatrix_1_2
    vardict = Dict{Int,Variable}()
    @RV vardict[3*i+1] ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @test vardict[3*i+1].id == :vardict_4

    # @RV should work with '= syntax'
    g = FactorGraph()
    @RV x = constant(1.0) + constant(2.0)
    @test length(g.variables) == 3
    @test x.id == :x
    @RV [id=:my_y] y = x + constant(2.0)
    @test length(g.variables) == 5
    @test y.id == :my_y

    # @RV without node definition should create a new Variable
    g = FactorGraph()
    @RV x
    @test isa(x, Variable)
    @test g.variables[:x] === x
    @RV [id=:x_new] x
    @test g.variables[:x_new] === x
    @test length(g.variables) == 2

    # @RV without node definition should create a new Variable with a given index
    g = FactorGraph()
    xi = Vector{Variable}(undef, 1)
    t = 1
    @RV xi[t]
    @test isa(xi[1], Variable)
    @test g.variables[:xi_1] === xi[1]

    # @RV should throw an exception on incorrect usage
    @test_throws ErrorException @RV 1
    @test_throws ErrorException @RV [ 1 ]
    @test_throws ErrorException @RV [ 1 ] 1
    @test_throws ErrorException @RV [ 1 ] x
    @test_throws ErrorException @RV [ x ] x
    @test_throws ErrorException @RV [ x = :x_id ]
    @test_throws ErrorException @RV [ x = :x_id ] 1
end

end #module
