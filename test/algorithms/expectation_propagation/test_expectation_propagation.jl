module ExpectationPropagationTest

using Base.Test
using ForneyLab
import ForneyLab: SoftFactor, generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable, EPSigmoidGP1, SPGaussianMeanVariancePPV

# Integration helper
type MockNode <: SoftFactor
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

@expectationPropagationRule(:node_type     => MockNode,
                            :outbound_type => Message{Gaussian},
                            :inbound_types => (Message{Gaussian}, Message{PointMass}),
                            :outbound_id   => 1,
                            :name          => EPMockGP1)

@testset "@expectationPropagationRule" begin
    @test EPMockGP1 <: ExpectationPropagationRule{MockNode}
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    m ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    nd = MockNode([m, constant(0.0)])
    inferred_outbound_types = Dict(nd.i[1].partner => Message{Gaussian}, nd.i[2].partner => Message{PointMass})

    entry = ScheduleEntry(nd.i[1], ExpectationPropagationRule{MockNode})
    inferUpdateRule!(entry, entry.msg_update_rule, inferred_outbound_types)

    @test entry.msg_update_rule == EPMockGP1
end

@testset "expectationPropagationSchedule" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    z = Variable[]
    nd_z = FactorNode[]
    for i = 1:3
        z_i = Variable()
        nd_z_i = Sigmoid(z_i, m)
        placeholder(z_i, :y, index=i)
        push!(z, z_i)
        push!(nd_z, nd_z_i)
    end

    schedule = expectationPropagationSchedule(m)

    @test length(schedule) == 15
    @test schedule[2] == ScheduleEntry(nd_z[1].i[:real], EPSigmoidGP1)
    @test schedule[4] == ScheduleEntry(nd_z[2].i[:real], EPSigmoidGP1)
    @test schedule[6] == ScheduleEntry(nd_z[3].i[:real], EPSigmoidGP1)
    @test schedule[10] == ScheduleEntry(nd_m.i[:out], SPGaussianMeanVariancePPV)
end

end # module