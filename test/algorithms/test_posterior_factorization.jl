module PosteriorFactorizationTest

using Test
using ForneyLab

import ForneyLab:Cluster

@testset "PosteriorFactorization" begin
    pfz = PosteriorFactorization()
    @test pfz.posterior_factors == Dict{Symbol, PosteriorFactor}()
    @test pfz.edge_to_posterior_factor == Dict{Edge, PosteriorFactor}()
    @test pfz.node_edge_to_cluster == Dict{Tuple{FactorNode, Edge}, Cluster}()
    @test currentPosteriorFactorization() === pfz
end

end # module