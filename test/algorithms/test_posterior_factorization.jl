module PosteriorFactorizationTest

using Test
using ForneyLab

using ForneyLab: Cluster, deterministicEdges, deterministicEdgeSet, isDeterministic

@testset "PosteriorFactorization" begin
    pfz = PosteriorFactorization()
    @test pfz.posterior_factors == Dict{Symbol, PosteriorFactor}()
    @test pfz.edge_to_posterior_factor == Dict{Edge, PosteriorFactor}()
    @test pfz.node_edge_to_cluster == Dict{Tuple{FactorNode, Edge}, Cluster}()
    @test currentPosteriorFactorization() === pfz
end

@testset "isDeterministic" begin
    # Clamp
    fg = FactorGraph()
    @RV x = constant(0.0)
    @test isDeterministic(fg.nodes[:clamp_1].i[:out], Dict{Interface, Bool}())
    
    # Stochastic
    fg = FactorGraph()
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    @test !isDeterministic(fg.nodes[:gaussianmeanvariance_1].i[:out], Dict{Interface, Bool}())

    # Equality
    fg = FactorGraph()
    @RV x = equal(constant(0.0), constant(0.0))
    @test isDeterministic(fg.nodes[:equality_1].i[1],
                          Dict{Interface, Bool}(fg.nodes[:clamp_1].i[:out] => true, 
                                                fg.nodes[:clamp_2].i[:out] => false))
    @test !isDeterministic(fg.nodes[:equality_1].i[1],
                           Dict{Interface, Bool}(fg.nodes[:clamp_1].i[:out] => false, 
                                                 fg.nodes[:clamp_2].i[:out] => false))
    # Deterministic
    fg = FactorGraph()
    @RV x ~ Addition(0.0, 0.0)
    @test isDeterministic(fg.nodes[:addition_1].i[:out],
                          Dict{Interface, Bool}(fg.nodes[:clamp_1].i[:out] => true, 
                                                fg.nodes[:clamp_2].i[:out] => true))
    @test !isDeterministic(fg.nodes[:addition_1].i[:out],
                           Dict{Interface, Bool}(fg.nodes[:clamp_1].i[:out] => true, 
                                                 fg.nodes[:clamp_2].i[:out] => false))
end

@testset "deterministicEdgeSet and deterministicEdges" begin
    fg = FactorGraph()
    
    @RV x ~ GaussianMeanVariance(0.0, 1.0)
    @RV y ~ GaussianMeanVariance(0.0, 1.0)
    @RV w = x + y
    @RV a = constant(1.0)
    @RV z = a*w
    placeholder(z, :z)

    e_x = x.edges[1]
    e_y = y.edges[1]
    e_w = w.edges[1]
    e_a = a.edges[1]    
    e_z = z.edges[1]
    e_1 = fg.nodes[:clamp_1].i[:out].edge
    e_2 = fg.nodes[:clamp_2].i[:out].edge
    e_3 = fg.nodes[:clamp_3].i[:out].edge
    e_4 = fg.nodes[:clamp_4].i[:out].edge
    
    @test deterministicEdgeSet(e_x) == Set{Edge}([e_w, e_a, e_z])
    @test deterministicEdges(fg) == Set{Edge}([e_w, e_a, e_z, e_1, e_2, e_3, e_4])
end

@testset "deterministicEdges for graph with dangling edge" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(1.0, 1.0)
    @RV y = 4.0 * x

    e_1 = fg.nodes[:clamp_1].i[:out].edge
    e_2 = fg.nodes[:clamp_2].i[:out].edge
    e_3 = fg.nodes[:clamp_3].i[:out].edge

    @test deterministicEdges(fg) == Set{Edge}([e_1, e_2, e_3])
end

end # module