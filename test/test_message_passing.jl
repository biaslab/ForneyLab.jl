module MessagePassingTest

using Test
using ForneyLab
using ForneyLab: generateId, addNode!, associate!, summaryPropagationSchedule, matches, flatten, Cluster

@testset "Message" begin
    msg = Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test isa(msg, Message{GaussianMeanVariance})
    @test isa(msg, Message{GaussianMeanVariance, Univariate})
    @test !isa(msg, Message{PointMass})
    @test !isa(msg, Message{GaussianMeanPrecision})
    @test !isa(msg, Message{GaussianMeanVariance, Multivariate})
    @test !isa(msg, Message{GaussianMeanVariance, MatrixVariate})
    @test msg.dist == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)

    @test Message(Univariate, PointMass, m=0.0) == Message(Univariate, PointMass, m=0.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, GaussianMeanVariance, m=0.0, v=1.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, PointMass, m=1.0)
end

@testset "matches" begin
    @test matches(Message{Gaussian}, Message)
    @test matches(Message{Gaussian, Univariate}, Message)
    @test matches(Message{Gaussian, Univariate}, Message{Gaussian, Univariate})
    @test matches(Message{GaussianMeanVariance, Univariate}, Message{Gaussian, Univariate})
    @test !matches(Message{GaussianMeanVariance, Univariate}, Message{Gaussian, Multivariate})
    @test !matches(Message{GaussianMeanVariance, Univariate}, Message{PointMass, Univariate})
    @test matches(Message{GaussianMeanVariance, Univariate}, Message{Gaussian})
    @test matches(Message{Gaussian}, Message{Gaussian})
    @test matches(Message{GaussianMeanVariance}, Message{Gaussian})
    @test !matches(Nothing, Message{Gaussian})
    @test !matches(Nothing, Message{GaussianMeanVariance})
    @test matches(Message{Gamma, Univariate}, Message{Union{Gamma, Wishart}, Univariate})
    @test matches(Message{Gamma}, Message{Union{Gamma, Wishart}})
    @test !matches(Message{GaussianMeanVariance, Univariate}, ProbabilityDistribution{Univariate, GaussianMeanVariance})
    @test !matches(ProbabilityDistribution{Univariate, GaussianMeanVariance}, Message{GaussianMeanVariance, Univariate})
end

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

@testset "ScheduleEntry" begin
    FactorGraph()
    nd = MockNode([Variable(), Variable()])

    # Schedule should be constructable by hand
    schedule = [ScheduleEntry(nd.i[1], Nothing), ScheduleEntry(nd.i[2], Nothing)]
    @test isa(schedule, Schedule)
end

@testset "summaryPropagationSchedule" begin
    # Construct
    #
    #    a    b    c
    # [1]--[2]--[3]--[4]
    #            |d
    #           [5]

    FactorGraph()

    a = Variable()
    b = Variable()
    c = Variable()
    d = Variable()

    n1 = MockNode([a])
    n2 = MockNode([a, b])
    n3 = MockNode([b, c, d])
    n4 = MockNode([c])
    n5 = MockNode([d])

    cl = Cluster(n1, [a.edges; b.edges])

    schedule = summaryPropagationSchedule(Variable[d], Cluster[cl])

    @test length(schedule) == 7
    @test ScheduleEntry(n1.i[1], Nothing) in schedule
    @test ScheduleEntry(n2.i[2], Nothing) in schedule
    @test ScheduleEntry(n4.i[1], Nothing) in schedule
    @test ScheduleEntry(n3.i[3], Nothing) in schedule
    @test ScheduleEntry(n5.i[1], Nothing) in schedule
    @test ScheduleEntry(n3.i[1], Nothing) in schedule
    @test ScheduleEntry(n2.i[1], Nothing) in schedule
end

@testset "flatten" begin
    FactorGraph()

    a = Variable()
    b = Variable()
    c = Variable()

    n1 = MockNode([a])
    n2 = MockNode([a, b])
    n3 = MockNode([b, c])
    n4 = MockNode([Variable()])

    schedule = [ScheduleEntry(n1.i[1], Nothing),
                ScheduleEntry(n2.i[2], Nothing),
                ScheduleEntry(n3.i[2], Nothing)]
    schedule[2].internal_schedule = [   ScheduleEntry(n1.i[1], Nothing),
                                        ScheduleEntry(n4.i[1], Nothing)]

    flat_schedule = flatten(schedule)

    @test length(flat_schedule) == 3
    @test flat_schedule[1] == schedule[1]
    @test flat_schedule[2] == schedule[2].internal_schedule[2]
    @test flat_schedule[3] == schedule[3]
end


end # module