module SumProductTest

using Base.Test
using ForneyLab
import ForneyLab: generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable
import ForneyLab: SPClamp

# Integration helper
type MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Array(Interface, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = associate!(Interface(self), vars[idx])
        end

        return self
    end
end

@sumProductRule(:node_type     => MockNode,
                :outbound_type => Message{Univariate{PointMass}},
                :inbound_types => (Void, Message{Univariate{PointMass}}, Message{Univariate{PointMass}}),
                :name          => SPMockOutPP)

@testset "@SumProductRule" begin
    @test SPMockOutPP <: SumProductRule{MockNode}
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    nd = MockNode([Variable(), constant(0.0), constant(0.0)])
    inferred_outbound_types = Dict(nd.i[2].partner => Message{Univariate{PointMass}}, nd.i[3].partner => Message{Univariate{PointMass}})

    entry = ScheduleEntry(nd.i[1], SumProductRule{MockNode})
    inferUpdateRule!(entry, entry.msg_update_rule, inferred_outbound_types)

    @test entry.msg_update_rule == SPMockOutPP
end

@testset "sumProductSchedule" begin
    FactorGraph()
    x = Variable()
    nd = MockNode([x, constant(0.0), constant(0.0)])

    schedule = sumProductSchedule(x)

    @test length(schedule) == 3
    @test ScheduleEntry(nd.i[2].partner, SPClamp) in schedule
    @test ScheduleEntry(nd.i[3].partner, SPClamp) in schedule
    @test ScheduleEntry(nd.i[1], SPMockOutPP) in schedule
end

end # module