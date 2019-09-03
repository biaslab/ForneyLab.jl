module JointMarginalsTest

using Test
using ForneyLab

import ForneyLab: generateId, addNode!, associate!, inferMarginalRule, isApplicable, Cluster, Product

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Vector{Interface}(undef, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = associate!(Interface(self), vars[idx])
        end

        return self
    end
end

@marginalRule(:node_type     => MockNode,
              :inbound_types => (Message{PointMass}, Message{PointMass}, ProbabilityDistribution),
              :name          => MMockPPD)

@testset "@marginalRule" begin
    @test MMockPPD <: MarginalRule{MockNode}
end

@testset "inferMarginalRule" begin
    FactorGraph()
    nd = MockNode([constant(0.0), constant(0.0), constant(0.0)])
    cluster = Cluster(nd, [nd.i[1].edge, nd.i[2].edge])

    @test inferMarginalRule(cluster, [Message{PointMass}, Message{PointMass}, ProbabilityDistribution]) == MMockPPD
end

@structuredVariationalRule(:node_type     => MockNode,
                           :outbound_type => Message{PointMass},
                           :inbound_types => (Nothing, Message{PointMass}, ProbabilityDistribution),
                           :name          => SVBMock1VGD)

@structuredVariationalRule(:node_type     => MockNode,
                           :outbound_type => Message{PointMass},
                           :inbound_types => (Message{PointMass}, Nothing, ProbabilityDistribution),
                           :name          => SVBMock2GVD)

@testset "marginalSchedule" begin
    FactorGraph()
    v1 = constant(0.0)
    v2 = constant(0.0)
    v3 = constant(0.0)
    nd = MockNode([v1, v2, v3])

    RecognitionFactorization()
    rf_12 = RecognitionFactor([v1, v2])

    schedule = variationalSchedule(rf_12)
    marginal_schedule = marginalSchedule(rf_12, schedule)

    @test length(marginal_schedule) == 1
    @test marginal_schedule[1].target == first(rf_12.clusters)
    @test marginal_schedule[1].interfaces[1] == nd.i[1].partner
    @test marginal_schedule[1].interfaces[2] == nd.i[2].partner
    @test marginal_schedule[1].marginal_update_rule == MMockPPD
end

end # module