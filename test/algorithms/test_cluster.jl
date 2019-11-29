module ClusterTest

using Test
using ForneyLab

import ForneyLab: Cluster

@testset "Cluster" begin
    g = FactorGraph()

    m = Variable(id=:m)
    v = Variable(id=:v)
    y = Variable(id=:y)
    nd = GaussianMeanVariance(y, m, v)
    em = nd.i[:m].edge
    ev = nd.i[:v].edge

    cluster = Cluster(nd, [em, ev])

    @test cluster.id == :m_v
    @test cluster.node == nd
    @test cluster.edges[1] == em
    @test cluster.edges[2] == ev
end

end # module
