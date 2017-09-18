module SumProductTest

using Base.Test
using ForneyLab
import ForneyLab: generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable, SPConstant

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
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Void),
                :name          => SPMockPPV)

@testset "@SumProductRule" begin
    @test SPMockPPV <: SumProductRule{MockNode}
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    nd = MockNode([constant(0.0), constant(0.0), Variable()])
    inferred_outbound_types = Dict(nd.i[1].partner => Message{PointMass}, nd.i[2].partner => Message{PointMass})

    entry = ScheduleEntry(nd.i[3], SumProductRule{MockNode})
    inferUpdateRule!(entry, entry.msg_update_rule, inferred_outbound_types)

    @test entry.msg_update_rule == SPMockPPV
end

@testset "sumProductSchedule" begin
    FactorGraph()
    x = Variable()
    nd = MockNode([constant(0.0), constant(0.0), x])

    schedule = sumProductSchedule(x)

    @test length(schedule) == 3
    @test schedule[1] == ScheduleEntry(nd.i[1].partner, SPConstant)
    @test schedule[2] == ScheduleEntry(nd.i[2].partner, SPConstant)
    @test schedule[3] == ScheduleEntry(nd.i[3], SPMockPPV)
end

end # module