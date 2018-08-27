module StructuredVariationalBayesTest

using Test
using ForneyLab
import ForneyLab: SoftFactor, generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable
import ForneyLab: VBGaussianMeanVarianceOut, SVBGaussianMeanPrecisionMGVD, SVBGaussianMeanPrecisionOutVGD, VBGammaOut, SVBGaussianMeanPrecisionW

# Integration helper
mutable struct MockNode <: SoftFactor
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

@structuredVariationalRule(:node_type     => MockNode,
                           :outbound_type => Message{PointMass},
                           :inbound_types => (Nothing, Message{PointMass}, ProbabilityDistribution),
                           :name          => SVBMock1VGD)

@structuredVariationalRule(:node_type     => MockNode,
                           :outbound_type => Message{PointMass},
                           :inbound_types => (Message{PointMass}, Nothing, ProbabilityDistribution),
                           :name          => SVBMock2GVD)

@structuredVariationalRule(:node_type     => MockNode,
                           :outbound_type => Message{PointMass},
                           :inbound_types => (ProbabilityDistribution, Nothing),
                           :name          => SVBMock3DV)

@testset "@structuredVariationalRule" begin
    @test SVBMock1VGD <: StructuredVariationalRule{MockNode}
    @test SVBMock2GVD <: StructuredVariationalRule{MockNode}
    @test SVBMock3DV <: StructuredVariationalRule{MockNode}
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    v1 = constant(0.0)
    v2 = constant(0.0)
    v3 = constant(0.0)
    nd = MockNode([v1, v2, v3])

    RecognitionFactorization()
    RecognitionFactor([v1, v2])
    RecognitionFactor(v3)

    entry1 = ScheduleEntry(nd.i[1], StructuredVariationalRule{MockNode})
    inferUpdateRule!(entry1, entry1.msg_update_rule, Dict{Interface, Type}(nd.i[2].partner => Message{PointMass}))
    @test entry1.msg_update_rule == SVBMock1VGD

    entry2 = ScheduleEntry(nd.i[2], StructuredVariationalRule{MockNode})
    inferUpdateRule!(entry2, entry2.msg_update_rule, Dict{Interface, Type}(nd.i[1].partner => Message{PointMass}))
    @test entry2.msg_update_rule == SVBMock2GVD

    entry3 = ScheduleEntry(nd.i[3], StructuredVariationalRule{MockNode})
    inferUpdateRule!(entry3, entry3.msg_update_rule, Dict{Interface, Type}())
    @test entry3.msg_update_rule == SVBMock3DV
end

@testset "variationalSchedule" begin
    g = FactorGraph()
    s_0 = Variable()
    nd_s_0 = GaussianMeanVariance(s_0, constant(0.0), constant(1.0))
    w = Variable()
    nd_w = Gamma(w, constant(1.0), constant(1.0))

    s = Variable[]
    nd_s = FactorNode[]

    s_min = s_0
    for i = 1:3
        s_i = Variable()
        nd_s_i = GaussianMeanPrecision(s_i, s_min, w)
        push!(s, s_i)
        push!(nd_s, nd_s_i)

        s_min = s_i
    end
    nd_s_i = GaussianMeanVariance(s_min, constant(0.0), constant(huge))
    push!(nd_s, nd_s_i)

    rf = RecognitionFactorization()
    q_w = RecognitionFactor(w)
    q_s = RecognitionFactor([s_0; s])

    schedule_q_s = variationalSchedule(q_s)
    @test length(schedule_q_s) == 8
    @test ScheduleEntry(nd_s_0.i[:out], VBGaussianMeanVarianceOut) in schedule_q_s
    @test ScheduleEntry(nd_s[4].i[:out], VBGaussianMeanVarianceOut) in schedule_q_s
    @test ScheduleEntry(nd_s[3].i[:m], SVBGaussianMeanPrecisionMGVD) in schedule_q_s
    @test ScheduleEntry(nd_s[2].i[:m], SVBGaussianMeanPrecisionMGVD) in schedule_q_s
    @test ScheduleEntry(nd_s[1].i[:m], SVBGaussianMeanPrecisionMGVD) in schedule_q_s
    @test ScheduleEntry(nd_s[1].i[:out], SVBGaussianMeanPrecisionOutVGD) in schedule_q_s
    @test ScheduleEntry(nd_s[2].i[:out], SVBGaussianMeanPrecisionOutVGD) in schedule_q_s
    @test ScheduleEntry(nd_s[3].i[:out], SVBGaussianMeanPrecisionOutVGD) in schedule_q_s

    schedule_q_w = variationalSchedule(q_w)
    @test length(schedule_q_w) == 6
    @test ScheduleEntry(nd_w.i[:out], VBGammaOut) in schedule_q_w
    @test ScheduleEntry(nd_s[1].i[:w], SVBGaussianMeanPrecisionW) in schedule_q_w
    @test ScheduleEntry(nd_s[2].i[:w], SVBGaussianMeanPrecisionW) in schedule_q_w
    @test ScheduleEntry(nd_s[3].i[:w], SVBGaussianMeanPrecisionW) in schedule_q_w
end

end # module