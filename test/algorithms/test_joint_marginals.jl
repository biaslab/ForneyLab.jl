module JointMarginalsTest

using Test
using ForneyLab
using ForneyLab: generateId, addNode!, associate!, inferMarginalRule, isApplicable, Cluster, Product, setTargets!, outboundType, variationalSchedule

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
              :inbound_types => (Message, Message, ProbabilityDistribution),
              :name          => MMockMMD)

@testset "@marginalRule" begin
    @test MMockMMD <: MarginalRule{MockNode}
    @test isApplicable(MMockMMD, [Message, Message, ProbabilityDistribution])
    @test !isApplicable(MMockMMD, [Message, Message, ProbabilityDistribution, ProbabilityDistribution])    
end

@testset "inferMarginalRule" begin
    FactorGraph()
    nd = MockNode([constant(0.0), constant(0.0), constant(0.0)])
    cluster = Cluster(nd, [nd.i[1].edge, nd.i[2].edge])

    @test inferMarginalRule(cluster, [Message, Message, ProbabilityDistribution]) == MMockMMD
end

mutable struct SPMockNode <: SumProductRule{MockNode} end
ForneyLab.outboundType(SPMockNode) = Message
ForneyLab.isApplicable(SPMockNode, input_types::Vector{<:Type}) = true

@testset "marginalTable" begin
    FactorGraph()
    v1 = Variable()
    v2 = Variable()
    v3 = constant(0.0)
    nd1 = MockNode([v1])
    nd2 = MockNode([v2])
    nd3 = MockNode([v1, v2, v3])

    pfz = PosteriorFactorization()
    pf_12 = PosteriorFactor(pfz)
    setTargets!(pf_12, pfz, external_targets=true)

    pf_12.schedule = variationalSchedule(pf_12)
    marginal_table = marginalTable(pf_12)

    @test length(marginal_table) == 1
    @test marginal_table[1].target == first(pf_12.target_clusters)
    @test marginal_table[1].interfaces[1] == nd1.i[1]
    @test marginal_table[1].interfaces[2] == nd2.i[1]
    @test marginal_table[1].marginal_update_rule == MMockMMD
end

end # module