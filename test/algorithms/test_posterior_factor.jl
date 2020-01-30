module PosteriorFactorTest

using Test
using ForneyLab

import ForneyLab: nodesConnectedToExternalEdges, Cluster, hasCollider, condense, flatten

@testset "PosteriorFactor" begin
    g = FactorGraph()
    @RV m ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV w ~ Gamma(constant(1.0), constant(1.0))
    y = Variable[]
    for i = 1:3
        @RV y_i ~ GaussianMeanPrecision(m, w)
        placeholder(y_i, :y, index=i)
        push!(y, y_i)
    end

    pfz = PosteriorFactorization()
    q_m = PosteriorFactor(m)
    @test q_m.id == :posteriorfactor_1
    @test q_m.variables == Set([m])
    @test q_m.clusters == Set{Cluster}()
    @test q_m.internal_edges == edges(m)
    @test pfz.posterior_factors[:posteriorfactor_1] === q_m

    q_w = PosteriorFactor(w)
    @test q_w.id == :posteriorfactor_2
    @test q_w.variables == Set([w])
    @test q_w.clusters == Set{Cluster}()
    @test q_w.internal_edges == edges(w)
    @test pfz.posterior_factors[:posteriorfactor_2] === q_w

    # Joint factorizations
    q_m_w = PosteriorFactor([m, w])
    @test q_m_w.id == :posteriorfactor_3
    @test q_m_w.variables == Set([m, w])
    @test length(q_m_w.clusters) == 3 
    @test q_m_w.internal_edges == edges(Set([m, w]))
    @test pfz.posterior_factors[:posteriorfactor_3] === q_m_w

    q_y = PosteriorFactor(y)
    @test q_y.id == :posteriorfactor_4
    @test q_y.variables == Set(y)
    @test q_y.clusters == Set{Cluster}()
    @test q_y.internal_edges == edges(Set(y))
    @test pfz.posterior_factors[:posteriorfactor_4] === q_y
end

@testset "hasCollider()" begin

    # [N]--->[+]<---[N]
    #     a   |   b
    #         v c
    #        [N]
    #         |
    #         v d
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = a + b
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    algo = InferenceAlgorithm()
    q = PosteriorFactor([a,b,c])
    @test hasCollider(q) == true

    # [N]--->[=]<---[N]
    #     a   |   b
    #         v c
    #        [N]
    #         |
    #         v d
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = equal(a, b)
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    algo = InferenceAlgorithm()
    q = PosteriorFactor([a,b,c])
    @test hasCollider(q) == true

    # [N]--->[=]--->[N]--->■
    #         |  a      b
    #         v 
    #        [N]
    #         |
    #         v c
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(a, 1.0)
    placeholder(b, :b)
    placeholder(c, :c)
    algo = InferenceAlgorithm()
    q = PosteriorFactor(a)
    @test hasCollider(q) == false

    # [N]--->[N]--->[N]--->■
    #     a      b      c

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(b, 1.0)
    placeholder(c, :c)
    algo = InferenceAlgorithm()
    q = PosteriorFactor([a,b])
    @test hasCollider(q) == false

    # [N]--->[N]--->■
    #     a      b 
    #
    # [N]--->[N]--->■
    #     c      d    

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(0.0, 1.0)
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(b, :b)
    placeholder(d, :d)
    algo = InferenceAlgorithm()
    q = PosteriorFactor([a,c])
    @test hasCollider(q) == false

    # [N]--->[+]<---[N]
    #     a   |   b
    #         v c
    # [N]<---[=]--->[N]
    #  |             |
    #  v d           v e 
    #  ■             ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = a + b
    @RV d ~ GaussianMeanVariance(c, 1.0)
    @RV e ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    placeholder(e, :e)
    algo = InferenceAlgorithm()
    q = PosteriorFactor([a,b,c])
    @test hasCollider(q) == true
end

end # module