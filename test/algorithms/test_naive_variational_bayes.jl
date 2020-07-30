module NaiveVariationalBayesTest

using Test
using ForneyLab
using ForneyLab: SoftFactor, generateId, addNode!, associate!, inferUpdateRule!, outboundType, isApplicable, setTargets!, messagePassingSchedule
using ForneyLab: VBGaussianMeanVarianceOut, VBGaussianMeanPrecisionM, SPEqualityGaussian

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

@naiveVariationalRule(:node_type     => MockNode,
                      :outbound_type => Message{PointMass},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBMockOut)

@testset "@naiveVariationalRule" begin
    @test VBMockOut <: NaiveVariationalRule{MockNode}
    @test isApplicable(VBMockOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution])
    @test !isApplicable(VBMockOut, [Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution])    
end

@testset "inferUpdateRule!" begin
    FactorGraph()
    nd = MockNode([Variable(), constant(0.0), constant(0.0)])

    entry = ScheduleEntry(nd.i[1], NaiveVariationalRule{MockNode})
    inferUpdateRule!(entry, entry.message_update_rule, Dict{Interface, Type}())

    @test entry.message_update_rule == VBMockOut
end

@testset "messagePassingSchedule" begin
    g = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    w = Variable()
    nd_w = Gamma(w, constant(1.0), constant(1.0))
    y = Variable[]
    nd_y = FactorNode[]
    for i = 1:3
        y_i = Variable()
        nd_y_i = GaussianMeanPrecision(y_i, m, w)
        placeholder(y_i, :y, index=i)
        push!(y, y_i)
        push!(nd_y, nd_y_i)
    end

    pfz = PosteriorFactorization()
    q_m = PosteriorFactor(m)

    setTargets!(q_m, pfz, external_targets=true)
    schedule = messagePassingSchedule(q_m)

    @test length(schedule) == 6
    @test ScheduleEntry(nd_m.i[:out], VBGaussianMeanVarianceOut) in schedule
    @test ScheduleEntry(nd_y[2].i[:m], VBGaussianMeanPrecisionM) in schedule
    @test ScheduleEntry(nd_y[3].i[:m], VBGaussianMeanPrecisionM) in schedule
    @test ScheduleEntry(nd_m.i[:out].partner.node.i[3].partner, SPEqualityGaussian) in schedule
    @test ScheduleEntry(nd_y[1].i[:m], VBGaussianMeanPrecisionM) in schedule
    @test ScheduleEntry(nd_m.i[:out].partner, SPEqualityGaussian) in schedule
end

@testset "messagePassingAlgorithm" begin
    fg = FactorGraph()
    m = Variable()
    nd_m = GaussianMeanVariance(m, constant(0.0), constant(1.0))
    w = Variable()
    nd_w = Gamma(w, constant(1.0), constant(1.0))
    y = Variable[]
    nd_y = FactorNode[]
    for i = 1:3
        y_i = Variable()
        nd_y_i = GaussianMeanPrecision(y_i, m, w)
        placeholder(y_i, :y, index=i)
        push!(y, y_i)
        push!(nd_y, nd_y_i)
    end

    pfz = PosteriorFactorization(m, w)
    algo = messagePassingAlgorithm(m)

    @test isa(algo, InferenceAlgorithm)
end

@testset "messagePassingAlgorithm with custom target" begin
    fg = FactorGraph()

    # Define model
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y = 2.0*x
    @RV z = y + 1.0
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    # Derive algorithm
    pfz = PosteriorFactorization(y)
    algo = messagePassingAlgorithm(y)
    code = algorithmSourceCode(algo)

    @test occursin("marginals[:y]", code)
end

end # module