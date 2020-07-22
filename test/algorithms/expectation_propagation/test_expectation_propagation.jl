module ExpectationPropagationTest

using Test
using ForneyLab
using ForneyLab: SoftFactor, generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable, setTargets!
using ForneyLab: EPProbitIn1GP, SPGaussianMeanVarianceOutNPP, SPClamp, VBGaussianMeanPrecisionOut, SPProbitOutNG, expectationPropagationSchedule, variationalExpectationPropagationSchedule

# Integration helper
mutable struct MockNode <: SoftFactor
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

@expectationPropagationRule(:node_type     => MockNode,
                            :outbound_type => Message{Gaussian},
                            :inbound_types => (Message{PointMass}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPMockIn1GP)

@testset "@expectationPropagationRule" begin
    @test EPMockIn1GP <: ExpectationPropagationRule{MockNode}
    @test isApplicable(EPMockIn1GP, [Message{PointMass}, Message{Gaussian}], 2)
    @test !isApplicable(EPMockIn1GP, [Message{PointMass}, Message{Gaussian}, Message{Gaussian}], 2)
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    @RV m ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    nd = MockNode([constant(0.0), m])
    inferred_outbound_types = Dict(nd.i[2].partner => Message{Gaussian}, nd.i[1].partner => Message{PointMass})

    entry = ScheduleEntry(nd.i[2], ExpectationPropagationRule{MockNode})
    inferUpdateRule!(entry, entry.message_update_rule, inferred_outbound_types)

    @test entry.message_update_rule == EPMockIn1GP
end

@testset "expectationPropagationSchedule" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    z = Variable[]
    nd_z = FactorNode[]
    for i = 1:3
        z_i = Variable()
        nd_z_i = Probit(z_i, m)
        placeholder(z_i, :y, index=i)
        push!(z, z_i)
        push!(nd_z, nd_z_i)
    end

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(pfz, target_variables=Set{Variable}([m]))
    setTargets!(pf, pfz)
    schedule = expectationPropagationSchedule(pf)

    @test length(schedule) == 15
    @test schedule[2] == ScheduleEntry(nd_z[2].i[:in1], EPProbitIn1GP)
    @test schedule[4] == ScheduleEntry(nd_z[3].i[:in1], EPProbitIn1GP)
    @test schedule[8] == ScheduleEntry(nd_m.i[:out], SPGaussianMeanVarianceOutNPP)
    @test schedule[11] == ScheduleEntry(nd_z[1].i[:in1], EPProbitIn1GP)
end

@testset "variationalExpectationPropagationSchedule" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    w = Variable()
    nd_w = Gamma(w, constant(0.01), constant(0.01))
    y = Variable()
    nd_y = GaussianMeanPrecision(y, m, w)
    z = Variable()
    nd_z = Probit(z, y)
    placeholder(z, :z)

    pfz = PosteriorFactorization()
    q_y_z = PosteriorFactor([y, z])
    setTargets!(q_y_z, pfz, external_targets=true)
    schedule = variationalExpectationPropagationSchedule(q_y_z)

    @test length(schedule) == 4
    @test ScheduleEntry(nd_y.i[:out], VBGaussianMeanPrecisionOut) in schedule
    @test ScheduleEntry(nd_z.i[:out].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(nd_z.i[:in1], EPProbitIn1GP) in schedule
    @test ScheduleEntry(nd_z.i[:out], SPProbitOutNG) in schedule
end

@testset "expectationPropagationAlgorithm" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    z = Variable[]
    nd_z = FactorNode[]
    for i = 1:3
        z_i = Variable()
        nd_z_i = Probit(z_i, m)
        placeholder(z_i, :y, index=i)
        push!(z, z_i)
        push!(nd_z, nd_z_i)
    end
    
    pfz = PosteriorFactorization()
    algo = expectationPropagationAlgorithm(m)

    @test isa(algo, InferenceAlgorithm)
end

@testset "variationalExpectationPropagationAlgorithm" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    w = Variable()
    nd_w = Gamma(w, constant(0.01), constant(0.01))
    y = Variable()
    nd_y = GaussianMeanPrecision(y, m, w)
    z = Variable()
    nd_z = Probit(z, y)
    placeholder(z, :z)

    pfz = PosteriorFactorization([y; z], m, w)
    algo = variationalExpectationPropagationAlgorithm(pfz)

    @test isa(algo, InferenceAlgorithm)    
end

end # module