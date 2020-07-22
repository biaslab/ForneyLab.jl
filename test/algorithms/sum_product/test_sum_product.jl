module SumProductTest

using Test
using ForneyLab
using ForneyLab: generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable, setTargets!, sumProductSchedule
using ForneyLab: SPClamp

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Array{Interface}(undef, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = associate!(Interface(self), vars[idx])
        end

        return self
    end
end

@sumProductRule(:node_type     => MockNode,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPMockOutPP)

@testset "@SumProductRule" begin
    @test SPMockOutPP <: SumProductRule{MockNode}
    @test isApplicable(SPMockOutPP, [Nothing, Message{PointMass}, Message{PointMass}])
    @test !isApplicable(SPMockOutPP, [Nothing, Message{PointMass}, Message{PointMass}, Message{PointMass}])    
end

# Composite definition for inferUpdateRule! testset
@composite TestComposite (b,a) begin
    @RV z ~ GaussianMeanVariance(a, constant(1.0))
    b = constant(2.0) + z
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    nd = MockNode([Variable(), constant(0.0), constant(0.0)])
    inferred_outbound_types = Dict(nd.i[2].partner => Message{PointMass}, nd.i[3].partner => Message{PointMass})

    entry = ScheduleEntry(nd.i[1], SumProductRule{MockNode})
    inferUpdateRule!(entry, entry.message_update_rule, inferred_outbound_types)

    @test entry.message_update_rule == SPMockOutPP

    # Internal msg passing tests
    FactorGraph()
    a = constant(0.0, id=:a)
    b = Variable(id=:b)
    tc = TestComposite(b, a)
    inferred_outbound_types = Dict(tc.i[:a].partner => Message{PointMass})
    entry = ScheduleEntry(tc.i[:b], SumProductRule{Nothing})

    inferUpdateRule!(entry, entry.message_update_rule, inferred_outbound_types)

    @test isdefined(entry, :internal_schedule)
    @test length(entry.internal_schedule) == 4
    @test entry.message_update_rule == entry.internal_schedule[end].message_update_rule
end

@testset "sumProductSchedule" begin
    FactorGraph()
    x = Variable()
    nd = MockNode([x, constant(0.0), constant(0.0)])

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(pfz, target_variables=Set{Variable}([x]))
    setTargets!(pf, pfz)
    schedule = sumProductSchedule(pf)

    @test length(schedule) == 3
    @test ScheduleEntry(nd.i[2].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(nd.i[3].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(nd.i[1], SPMockOutPP) in schedule
end

@testset "sumProductAlgorithm" begin
    FactorGraph()
    x = Variable()
    nd = MockNode([x, constant(0.0), constant(0.0)])

    PosteriorFactorization()
    algo = sumProductAlgorithm(x)

    @test isa(algo, InferenceAlgorithm)
end

end # module