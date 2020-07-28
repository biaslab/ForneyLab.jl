module PosteriorFactorTest

using Test
using ForneyLab

using ForneyLab: nodesConnectedToExternalEdges, Cluster, condense, flatten, setTargets!

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
    @test q_m.internal_edges == edges(m)
    @test pfz.posterior_factors[:posteriorfactor_1] === q_m

    q_w = PosteriorFactor(w)
    @test q_w.id == :posteriorfactor_2
    @test q_w.internal_edges == edges(w)
    @test pfz.posterior_factors[:posteriorfactor_2] === q_w

    # Joint factorizations
    q_m_w = PosteriorFactor([m, w])
    @test q_m_w.id == :posteriorfactor_3
    @test q_m_w.internal_edges == edges(Set([m, w]))
    @test pfz.posterior_factors[:posteriorfactor_3] === q_m_w

    q_y = PosteriorFactor(y)
    @test q_y.id == :posteriorfactor_4
    @test q_y.internal_edges == edges(Set(y))
    @test pfz.posterior_factors[:posteriorfactor_4] === q_y
end

@testset "setTargets!" begin
    fg = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(0.0, 1.0)
    @RV z = x + y
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([z]))
    @test pfz.free_energy_flag == false
    @test pf.target_variables == Set{Variable}([z])
    @test pf.target_clusters == Set{Cluster}()

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([z]), external_targets=true)
    @test pfz.free_energy_flag == false
    @test pf.target_variables == Set{Variable}([z])
    @test pf.target_clusters == Set{Cluster}()

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([z]), free_energy=true)
    @test pfz.free_energy_flag == true
    @test pf.target_variables == Set{Variable}([x, y, z])
    @test length(pf.target_clusters) == 1
    @test first(pf.target_clusters).id == :x_y
end

end # module